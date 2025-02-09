/*		Copyright 1994 by Daniel R. Grayson		*/

#include "expr-exports.h"

#include "M2mem.h"
#include "types.h"
#include "debug.h"
#include "supervisorinterface.h"

#include <readline/readline.h>
#include <readline/history.h>
#ifdef WITH_MPI
#include <mpi.h>
#endif
#include <stdio.h>

extern struct JumpCell abort_jmp;
extern struct JumpCell interrupt_jmp;

extern void alarm_handler(int sig);
extern void interrupt_handler(int sig);
extern void oursignal(int sig, void (*handler)(int));

int have_arg_no_int = 0;

// TODO: remove this from actors5.d and system.d?
void system_handleInterruptsSetup(M2_bool handleInterrupts) {
  if (!have_arg_no_int) {
    oursignal(SIGINT,handleInterrupts ? interrupt_handler : SIG_DFL);
  }
  oursignal(SIGALRM,handleInterrupts ? alarm_handler : SIG_DFL);
}

static double startTime;
double system_cpuTime(void) {
  struct timespec t;
  int err = clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t);
  if (err) return 0; /* silent about error */
  double u = t.tv_sec + t.tv_nsec * 1e-9;
  return u - startTime;
}
void system_cpuTime_init(void) {
  startTime = system_cpuTime();
}

// begin MPI ------------------------------------------------------------------
#ifdef WITH_MPI
int MPInumberOfProcesses() {
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  return world_size;
}

int MPImyProcessNumber() {
  int world_rank;
  if (MPI_Comm_rank(MPI_COMM_WORLD, &world_rank) != MPI_SUCCESS)
    world_rank = -1;
  return world_rank;
}

// blocking send 
int MPIsendString(M2_string s, int p) {
  char *t = M2_tocharstar(s);
  int ret = MPI_Send(t, strlen(t)+1, MPI_CHAR, p, 0 /*tag*/, MPI_COMM_WORLD);
  GC_FREE(t); 
  return ret;
}

// blocking receive
M2_string MPIreceiveString(int p) {
  MPI_Status status;
  // Probe for an incoming message from process zero
  MPI_Probe(p, 0/*tag*/, MPI_COMM_WORLD, &status);
  // When probe returns, the status object has the size and other
  // attributes of the incoming message. Get the message size
  int size;
  MPI_Get_count(&status, MPI_CHAR, &size);
  // Allocate a buffer to hold the incoming numbers
  char* s = (char*) malloc(sizeof(char) * size);
  // Now receive the message with the allocated buffer
  MPI_Recv(s, size, MPI_CHAR, p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  M2_string ret = M2_tostring(s);
  free(s);
  return ret;
}


// nonblocking probe
int MPIprobeInterrupt() {
  int flag = 0;
  int receive_tag;
  MPI_Status status;
  const int master = 0;
  //!!! consider any integer message an interrupt???
  return MPI_Iprobe(master, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
}

// nonblocking send (needed???)
int MPIsendStringNonblocking(M2_string s, int p) {
  char *t = M2_tocharstar(s);
  MPI_Request request;
  int ret = MPI_Isend(t, strlen(t)+1, MPI_CHAR, p, 0 /*tag*/, MPI_COMM_WORLD, &request);
  // should remember "request" if e.g. need to check the completion 
  GC_FREE(t); // is it OK to free? (Should be if MPI still refers to this... but does it?)  
  return ret;
}

// interrupt
int MPIinterrupt(int p) {
  const int MPI_INTERRUPT_TAG = 3210; //!!! is also in e/interrupted.cpp
  int ret = MPI_Send(NULL, 0, MPI_INT, p, MPI_INTERRUPT_TAG, MPI_COMM_WORLD);
  // should remember "request" if e.g. need to check the completion 
  return ret;
}

#endif
// end MPI --------------------------------------------------------------------

double system_threadTime(void) {
  struct timespec t;
  int err = clock_gettime(CLOCK_THREAD_CPUTIME_ID, &t);
  if (err) return 0; /* silent about error */
  return t.tv_sec + t.tv_nsec * 1e-9;
}

void clean_up(void) {
  extern void close_all_dbms();
  close_all_dbms();
  while (pre_final_list != NULL) {
    pre_final_list->func();
    pre_final_list = pre_final_list->next;
  }
  while (final_list != NULL) {
    final_list->func();
    final_list = final_list->next;
  }
#ifdef WITH_PYTHON
  if (Py_IsInitialized()) Py_Finalize();
#endif
#ifndef NDEBUG
  trap();
#endif
#ifdef WITH_MPI
  int n = MPImyProcessNumber();
  if (n>0)
      printf("MPI: Bye world from process %d out of %d processes\n",
	     MPImyProcessNumber(), MPInumberOfProcesses());
  MPI_Finalize();
#endif
}

int system_isReady(int fd) {
  int ret;
  static fd_set r, w, e;
  struct timeval timeout;
  FD_SET(fd,&r);
  timerclear(&timeout);
  ret = select(fd+1,&r,&w,&e,&timeout);
  FD_CLR(fd,&r);
  return ret;
}

int system_hasException(int fd) {
  int ret;
  static fd_set r, w, e;
  struct timeval timeout;
  FD_SET(fd,&e);
  timerclear(&timeout);
  ret = select(fd+1,&r,&w,&e,&timeout);
  FD_CLR(fd,&e);
  return ret;
}

extern void M2_stack_trace();

extern void fatal(const char *s, ...);

extern void fatalarrayindex(int indx, int len, const char *file, int line, int column);

extern void fatalarraylen(int len, const char *file, int line, int column);

/******************************************************************************/
/*  Functions dealing with libreadline and completion                         */
/******************************************************************************/

static char *M2_completion_generator(const char *text, int state) {
  static int i;
  static char **v;
  char *p;
  if (state == 0) {
    M2_string s;
    M2_ArrayString ret;
    i = 0;
#ifdef free
#warning "'free' defined as macro, but we want to use the libc function"
#define free x
#endif
    if (v != NULL) free(v);
    s = M2_tostring(text);
    ret = expr_completions(s);
    freemem(s);
    v = M2_tocharstarstarmalloc(ret); /* readline will use free() to free these strings */
    freemem(ret);
  }
  p = v[i];
  if (p != NULL) i++;
  return p;
}

static char **M2_completion(const char *text, int start, int end) {
  rl_attempted_completion_over = TRUE;
  /* if (start > 0 && rl_line_buffer[start-1] == '"') ... filename completion ... */
  return rl_completion_matches(text, M2_completion_generator);
}

int system_readHistory(char *filename) { return read_history(filename); }
int system_appendHistory(int n, char *filename)
{
  return append_history(n, filename);
}

void system_addHistory(char *buf) { add_history(buf); }
char *system_getHistory(const int n)
{
  HIST_ENTRY *entry = history_get(n);
  if (entry != NULL) return entry->line;
  return NULL;
}

int system_historyLength() { return history_length; }

void system_initReadlineVariables(void) {
  static char readline_name[] = "M2";
  static char basic_word_break_characters[] = "!\"#$%&'()*+,-./:;<=>?@[\\]^_`{|}~ \t\n\r";
  rl_catch_signals = FALSE; /* tell readline not to catch signals, such as SIGINT */
  rl_readline_name = readline_name;
  rl_attempted_completion_function = M2_completion;
  rl_basic_word_break_characters = basic_word_break_characters;
  using_history();		/* this might also initialize readine, by calling rl_readline, on Mac OS X */
}

static int read_via_readline(char *buf,int len,char *prompt) {
  static char *p;		/* buffer, NULL if newline has already been returned */
  static int plen;		/* number of chars in p */
  static int i;			/* number of chars in p already returned */
  int r;			/* number of chars to return this time */
  if (len == 0) return 0;
  if (p == NULL) {
    interrupt_jmp.is_set = TRUE; /* for the interrupt handler */
    if (SETJMP(interrupt_jmp.addr)) { /* long jump occurred */
	 fprintf(stderr,"^C\n");
	 interrupt_jmp.is_set = FALSE;
	 rl_cleanup_after_signal();
	 rl_free_line_state();
	 return ERROR;
	 }
    p = readline(prompt);
    interrupt_jmp.is_set = FALSE;
    if (p == NULL) return 0;	/* EOF */
    i = 0;
    plen = strlen(p);
    add_history(p);
  }
  r = plen - i;
  if (r > len) r = len;
  memmove(buf,p+i,r), i+=r;
  if (i == plen && r < len) {
    free(p), p = NULL;
    buf[r++] = '\n';		/* readline() doesn't include the \n at the end */
  }
  return r;
}

int system_readline(M2_string buffer, int len, int offset, M2_string prompt) {
  char *p = M2_tocharstar(prompt);
  int r;
  if (offset < 0 || (int)buffer->len - offset < len) fatalarrayindex(len,buffer->len,__FILE__,__LINE__,-1);
  r = read_via_readline((char *)buffer->array + offset,len,p);
  freemem(p);
  return r;
}

/*
// Local Variables:
// compile-command: "echo \"make: Entering directory \\`$M2BUILDDIR/Macaulay2/d'\" && make -C $M2BUILDDIR/Macaulay2/d M2lib.o "
// tags-file-name: "TAGS"
// End:
*/

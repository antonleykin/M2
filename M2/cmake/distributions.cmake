set(CMAKE_VERBOSE_MAKEFILE OFF)

set(PACKAGE_TARNAME ${PROJECT_NAME})
set(DISTRIBUTION 1) # use this starting number to sequentially number the downstream distributions

## Setting variables for installation directories
## This is using https://cmake.org/cmake/help/latest/module/GNUInstallDirs.html#module:GNUInstallDirs
## Which follows https://www.gnu.org/prep/standards/html_node/Directory-Variables.html
## Also see FHS  https://refspecs.linuxfoundation.org/FHS_3.0/fhs/index.html
include(GNUInstallDirs)

# FIXME: usr-dist/common/share gives error?!
set(SHARE_DIR ${CMAKE_BINARY_DIR}/Macaulay2/share)
set(DIST_DIR  ${CMAKE_BINARY_DIR}/usr-dist)
set(CORE_DIR  ${SHARE_DIR}/Macaulay2/Core)

################################################################
# TODO
## Here we define the Macaulay2 layout, once and for all.
## There is a hidden dependency: see M2/distributions/dmg/Makefile.in and adjust the number "../"s if the basic layout here is changed.

set(prefix	${CMAKE_INSTALL_PREFIX}/common)
set(exec_prefix	${CMAKE_INSTALL_PREFIX}/${MACHINE})

set(common	"\$\{prefix\}")
set(exec	"\$\{exec_prefix\}")

set(datarootdir	${prefix}/share)
set(docdir	${datarootdir}/doc/${PACKAGE_TARNAME})
set(htmldir	${docdir}) # TODO: use this for documentation?

set(bindir	${exec}/${CMAKE_INSTALL_BINDIR})
set(infodir	${common}/${CMAKE_INSTALL_INFODIR})
set(datadir	${common}/${CMAKE_INSTALL_DATADIR})
set(mandir	${common}/${CMAKE_INSTALL_MANDIR})
set(docdir	${common}/${CMAKE_INSTALL_DOCDIR})
set(libdir	${exec}/${CMAKE_INSTALL_LIBDIR})
set(libexecdir	${exec}/${CMAKE_INSTALL_LIBEXECDIR})

set(gftablesdir	${datadir}/Macaulay2/Core/factory/)
set(packagesdir	${datadir}/Macaulay2)
set(libm2dir	${libdir}/Macaulay2)
set(emacsdir	${datadir}/emacs/site-lisp)
set(librariesdir ${libm2dir}/lib)
set(programsdir	${libexecdir}/Macaulay2/bin)
set(licensesdir	${libexecdir}/Macaulay2/program-licenses)
set(packagecachecoredir ${libm2dir}/cache)

###################################################################################
### STILL PROCESSING

#set(tail_packagesdir	${tail_datadir}/Macaulay2)
#set(tail_libm2dir	${tail_libdir}/Macaulay2)
#set(tail_emacsdir	${tail_datadir}/emacs}/site-lisp)
#set(tail_librariesdir	${tail_libm2dir}/lib)
#set(tail_programsdir	${tail_libexecdir}/Macaulay2}/bin)
#set(tail_licensesdir	${tail_libexecdir}/Macaulay2}/program-licenses)
#set(tail_packagecachecoredir	${tail_libm2dir}/cache)

#set(pre_packagesdir	${pre_datadir}/Macaulay2)
#set(pre_libm2dir	${pre_libdir}/Macaulay2)
#set(pre_emacsdir	${pre_datadir}/emacs/site-lisp)
#set(pre_librariesdir	${pre_libm2dir}/lib)
#set(pre_programsdir	${pre_libexecdir}/Macaulay2/bin)
#set(pre_licensesdir	${pre_libexecdir}/Macaulay2/program-licenses)
#set(pre_packagecachecoredir	${pre_libm2dir}/cache)


## We put docdir and datarootdir at end, because we are changing the values of the variables, and
## other variables depend on them.

## We don't normalize sysconfdir, sharedstatedir, and localstatedir, because Fedora
## insists in /usr/share/config.site that they be /etc, /var, and /var, respectively.

## We don't use those directories, anyway.  To ensure that, we set them to /nowhere.
## (We could avoid config.site with the CONFIG_SITE environment variable.)
#sysconfdir=/nowhere
#sharedstatedir=/nowhere
#localstatedir=/nowhere

#for i in bindir datadir includedir infodir libdir libexecdir mandir sbindir \
#         psdir pdfdir dvidir htmldir localedir gftablesdir docdir datarootdir
#do eval w=\$$i ; eval v="$w" ; eval v="$v" ; eval v="$v" ; eval v="$v" ; eval v="$v"
#   case $v in
#     "$exec_prefix"|"$exec_prefix"/*)
#	   eval tail_${i}=`echo $v | sed s,"^$exec_prefix/","",`
#	   eval  pre_${i}=`echo $v | sed s,"^$exec_prefix","'\\${pre_exec_prefix}'",`
#	   eval      ${i}=`echo $v | sed s,"^$exec_prefix","'\\${exec_prefix}'",`
#	   ;;
#     "$prefix"|"$prefix"/*)
#	   eval tail_${i}=`echo $v | sed s,"^$prefix/","",`
#	   eval  pre_${i}=`echo $v | sed s,"^$prefix","'\\${pre_prefix}'",`
#	   eval      ${i}=`echo $v | sed s,"^$prefix","'\\${prefix}'",`
#	   ;;
#     *) if test $i != gftablesdir #only needs normalized if we build factory
#	then AC_MSG_ERROR([expected "\${$i}" => "$w" to start with "\${prefix}" or with "\${exec_prefix}"])
#	fi ;;
#   esac
#done
#prefix=$save_prefix
#exec_prefix=$save_exec_prefix

#AC_SUBST(pre_bindir) dnl In Macaulay2/bin/M2.in we assume that $pre_bindir/.. is $pre_exec_prefix, which follows from pre_bindir=${pre_exec_prefix}/bin.
#AC_SUBST(pre_datadir)
#AC_SUBST(pre_includedir)
#AC_SUBST(pre_infodir)
#AC_SUBST(pre_libdir)
#AC_SUBST(pre_libexecdir)
#AC_SUBST(pre_localstatedir)
#AC_SUBST(pre_mandir)
#AC_SUBST(pre_sbindir)
#AC_SUBST(pre_sharedstatedir)
#AC_SUBST(pre_sysconfdir)
#AC_SUBST(pre_psdir)
#AC_SUBST(pre_pdfdir)
#AC_SUBST(pre_dvidir)
#AC_SUBST(pre_htmldir)
#AC_SUBST(pre_localedir)
#AC_SUBST(pre_docdir)
#AC_SUBST(pre_datarootdir)

#AC_SUBST(tail_bindir)
#AC_SUBST(tail_datadir)
#AC_SUBST(tail_includedir)
#AC_SUBST(tail_infodir)
#AC_SUBST(tail_libdir)
#AC_SUBST(tail_libexecdir)
#AC_SUBST(tail_localstatedir)
#AC_SUBST(tail_mandir)
#AC_SUBST(tail_sbindir)
#AC_SUBST(tail_sharedstatedir)
#AC_SUBST(tail_sysconfdir)
#AC_SUBST(tail_psdir)
#AC_SUBST(tail_pdfdir)
#AC_SUBST(tail_dvidir)
#AC_SUBST(tail_htmldir)
#AC_SUBST(tail_localedir)
#AC_SUBST(tail_docdir)
#AC_SUBST(tail_datarootdir)

#AC_SUBST(COMMONSTAGINGAREA,usr-dist) # staging area for common files
#AC_SUBST(LOCALSTAGINGAREA,usr-dist) # staging area for arch. dep. files
#AC_ARG_WITH(staging-area,
#    AS_HELP_STRING(--with-staging-area=...,directory for pre-installation of architecture-independent files (usr-dist)),
#    COMMONSTAGINGAREA=$withval)
#AC_ARG_ENABLE(common-staging-area,
#    AS_HELP_STRING(--enable-common-staging-area,use the common staging area),
#    if test "$enableval" = yes
#    then if test "$COMMONSTAGINGAREA" != usr-dist
#	 then AC_MSG_ERROR(--with-staging-area and --enable-common-staging-area options both provided)
#	 fi
#	 COMMONSTAGINGAREA=\${abs_top_srcdir}/BUILD/CommonStagingArea # this depends on config.Makefile.in setting abs_top_srcdir
#    fi
#    )
#[
#case `pwd` in
#  .) ;;
#  *)
#    case "$COMMONSTAGINGAREA" in
#      [\\/]* | ?:[\\/]* | '${abs'* );;
#      .) COMMONSTAGINGAREA=`pwd`;;
#      *) COMMONSTAGINGAREA=`pwd`/"$COMMONSTAGINGAREA";;
#    esac
#    case "$LOCALSTAGINGAREA" in
#      [\\/]* | ?:[\\/]* | '${abs'* );;
#      .) LOCALSTAGINGAREA=`pwd`;;
#      *) LOCALSTAGINGAREA=`pwd`/"$LOCALSTAGINGAREA";;
#    esac;;
#esac
#]
#AC_SUBST(pre_prefix,$COMMONSTAGINGAREA/common) # as in layout.m2.in
#AC_SUBST(pre_exec_prefix,$LOCALSTAGINGAREA/$MACHINE) # as in layout.m2.in
#AC_MSG_NOTICE([staging area for common files: $pre_prefix])
#AC_MSG_NOTICE([staging area for architecture dependent files: $pre_exec_prefix])



################################################################

#AC_OUTPUT() # distribution files configured

#distributions/top/INSTALL
#distributions/top/preremove
#distributions/top/postinstall
#distributions/rpm/Macaulay2-body.spec
#distributions/rpm/Macaulay2-common-body.spec
#distributions/cygwin/Macaulay2.README
#distributions/cygwin/setup.hint-base
#distributions/cygwin/postinstall.sh
#distributions/cygwin/preremove.sh
#distributions/cygwin/server/Macaulay2-icons/setup.hint
#distributions/freebsd/post-install
#distributions/freebsd/post-deinstall
#distributions/deb/macaulay2-common/postrm
#distributions/deb/macaulay2-common/postinst
#distributions/deb/macaulay2-common/prerm
#distributions/deb/macaulay2-common/preinst
#distributions/deb/macaulay2/postrm
#distributions/deb/macaulay2/postinst
#distributions/deb/macaulay2/prerm
#distributions/deb/macaulay2/preinst

#AC_SUBST(COMPRESS,gz)  AC_ARG_ENABLE(compress, AS_HELP_STRING([--enable-compress=[gz|bz2]],compression method for tarball), COMPRESS=$enableval)

#AC_SUBST(M2TARFILE,no) AC_ARG_ENABLE(tarfile, AS_HELP_STRING(--enable-tarfile,prepare binary and source packages as compressed tar files), M2TARFILE=$enableval)
#AC_SUBST(TARLIBS,no)   AC_ARG_ENABLE(tarlibs, AS_HELP_STRING(--enable-tarlibs,include symbolic links to needed shared libraries for tar), TARLIBS=$enableval)

#AC_SUBST(DEB,no)
#AC_ARG_ENABLE(deb,  AS_HELP_STRING(--enable-deb,prepare a *.deb package (for debian, ubuntu, ...)), DEB=$enableval)
#if test "$DEB" = yes
#then if test "$OPTIMIZE" = no
#     then AC_MSG_ERROR([--disable-optimize and --enable-deb both specified])
#     fi
#     if test "$DEBUG" = yes
#     then AC_MSG_ERROR([--enable-debug and --enable-deb both specified])
#     fi
#fi

#AC_SUBST(FREEBSD,no)
#AC_ARG_ENABLE(freebsd,  AS_HELP_STRING(--enable-freebsd,prepare a package file for freebsd), FREEBSD=$enableval)
#if test "$FREEBSD" = yes
#then if test "$OPTIMIZE" = no
#     then AC_MSG_ERROR([--disable-optimize and --enable-freebsd both specified])
#     fi
#     if test "$DEBUG" = yes
#     then AC_MSG_ERROR([--enable-debug and --enable-freebsd both specified])
#     fi
#fi

#AC_SUBST(RPM,no)
#AC_ARG_ENABLE(rpm,  AS_HELP_STRING(--enable-rpm,prepare a *.rpm package (for red hat based systems)), RPM=$enableval)
#if test "$RPM" = yes
#then if test "$OPTIMIZE" = no
#     then AC_MSG_ERROR([--disable-optimize and --enable-rpm both specified])
#     fi
#     if test "$DEBUG" = yes
#     then AC_MSG_ERROR([--enable-debug and --enable-rpm both specified])
#     fi
#fi

#AC_SUBST(DMG,no)
#AC_ARG_ENABLE(dmg,  AS_HELP_STRING(--enable-dmg,prepare a *.dmg package (for Mac OS)), DMG=$enableval)
#if test "$DMG" = yes
#then if test "$OPTIMIZE" = no
#     then AC_MSG_ERROR([--disable-optimize and --enable-dmg both specified])
#     fi
#     if test "$DEBUG" = yes
#     then AC_MSG_ERROR([--enable-debug and --enable-dmg both specified])
#     fi
#fi

################################################################

AC_INIT(MAQ, 0.7.1)
AM_CONFIG_HEADER(config.h)
AM_INIT_AUTOMAKE([no-dependencies])

AC_CANONICAL_HOST
AC_PROG_CC
AC_PROG_CXX

# set CFLAGS and CXXFLAGS
user_CFLAGS=${CFLAGS}
generic_CFLAGS="-Wall"
ext_CFLAGS=""
case "${host_cpu}-${host_os}" in
  i386-darwin*)
    AC_ARG_ENABLE(universal, [  --enable-universal      enable universal binary (Intel Apple only)],
				  [is_uni=1], [is_uni=0])
	case $is_uni in
	  1) CFLAGS="-m64"
	     AC_COMPILE_IFELSE([AC_LANG_PROGRAM],
						   [ext_CFLAGS="-arch x86_64 -arch i386 -arch ppc64 -arch ppc"],
						   [ext_CFLAGS="-arch i386 -arch ppc"]);;
	  0) CFLAGS="-m64"
	  	 AC_COMPILE_IFELSE([AC_LANG_PROGRAM], [ext_CFLAGS="-m64"], []);;
	esac;;
  *)
    AC_MSG_CHECKING([if gcc accepts -m64])
    CFLAGS="-m64"
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM], [ext_CFLAGS="-m64"; AC_MSG_RESULT([yes])],
					  [ext_CFLAGS="-D_FILE_OFFSET_BITS=64"; AC_MSG_RESULT([no])]);;
esac
AC_ARG_ENABLE(experimental, [  --enable-experimental   enable experimental features],
			  [ext_CFLAGS="${ext_CFLAGS} -DMAQ_SHOW_EXPERIMENTAL"], [])
AC_ARG_ENABLE(slowmap,      [  --enable-slowmap        use slowmap mode for better sensitivity],
			  [], [ext_CFLAGS="${ext_CFLAGS} -D_FASTMAP"])
AC_ARG_ENABLE(shortread,    [  --enable-shortreads     use shortread mode],
			  [], [ext_CFLAGS="${ext_CFLAGS} -DMAQ_LONGREADS"])
AC_ARG_ENABLE(intel64,      [  --enable-intel64        optimize for Intel64 CPU such as Xeon and Core2],
			  [ext_CFLAGS="${ext_CFLAGS} -mtune=nocona"], [])
CFLAGS="${generic_CFLAGS} ${ext_CFLAGS} ${user_CFLAGS}"
CXXFLAGS=$CFLAGS

AC_STDC_HEADERS
AC_CHECK_HEADER(zlib.h)

AC_OUTPUT(Makefile)

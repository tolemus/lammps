
maise_INC="-I../MAISE/inc -I../MAISE/lib/include"
maise_LIB="-lmaise -lgomp -L../MAISE/lib"

if (test $1 = 1) then

  if (test -e ../Makefile.package) then
    sed -i -e "s|^PKG_INC.=[ \t]*|&$maise_INC |" ../Makefile.package
    sed -i -e "s|^PKG_LIB.=[ \t]*|&$maise_LIB |" ../Makefile.package
  fi
elif (test $1 = 0) then
  if (test -e ../Makefile.package) then
    sed -i -e "s|\(^PKG_INC.=.*\)\($maise_INC\)\(.*\)|\1\3|" ../Makefile.package
    sed -i -e "s|\(^PKG_LIB.=.*\)\($maise_LIB\)\(.*\)|\1\3|" ../Makefile.package
  fi
else
  echo Unknown \$1 $1
fi

maise_INC="-I../MAISE/inc -I../MAISE/lib/include"
maise_LIB="-lmaise -lgomp -L../MAISE/lib"

if (test $1 = 1) then

  if (test -e ../Makefile.package) then
    sed -i -e "s|^PKG_INC.=[ \t]*|&$maise_INC |" ../Makefile.package
    sed -i -e "s|^PKG_LIB.=[ \t]*|&$maise_LIB |" ../Makefile.package
  fi
else
    echo Unknown \$1 $1
fi
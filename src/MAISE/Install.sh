
maise_MAKE="include ../MAISE/Makefile.maise"

if (test $1 = 1) then
  # Install
  sed -i -e "\|$maise_MAKE|d" ../Makefile.package.settings
  sed -i -e 's|\(^PKG_INC.=.*\)\($(maise_INC)\)\(.*\)|\1\3|' ../Makefile.package
  sed -i -e 's|\(^PKG_LIB.=.*\)\($(maise_LIB)\)\(.*\)|\1\3|' ../Makefile.package
  cp ../MAISE/*_maise.* ..
  if (test -e ../Makefile.package) then
    sed -i -e "\|$maise_MAKE|d" ../Makefile.package.settings
    echo $maise_MAKE >> ../Makefile.package.settings
    sed -i -e 's|^PKG_INC.=[ \t]*|&$(maise_INC) |' ../Makefile.package
    sed -i -e 's|^PKG_LIB.=[ \t]*|&$(maise_LIB) |' ../Makefile.package
  fi
elif (test $1 = 0) then
  # Uninstall
  rm ../*_maise.*
  # if (test -e ../Makefile.package) then
  # fi
else
  echo Unknown \$1 $1
fi
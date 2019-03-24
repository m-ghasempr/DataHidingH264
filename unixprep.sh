#/bin/sh

echo "Creating obj directories..."
if test -d lencod/obj ;then echo "  lencod/obj exists." ; else mkdir lencod/obj ;fi
if test -d ldecod/obj ;then echo "  ldecod/obj exists." ; else mkdir ldecod/obj ;fi

if test -z "`which sed`" ;then
  echo "Fatal: \"sed\" not found!"
  exit
fi

touch lencod/dependencies
touch ldecod/dependencies

echo "Removing DOS LF chars..."
for s in `find . -name \*.[ch]`
do
cat $s|sed -e "s///" >$s.tmp
mv $s.tmp $s
done
for s in `find . -name Makefile`
do
cat $s|sed -e "s///" >$s.tmp
mv $s.tmp $s
done

echo "Done."


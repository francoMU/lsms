
mv k.out k.out.new
diff k.out.new k.out.ref > /dev/null 2>&1
error=$?
if [ $error -eq 0 ]
then
   exit 0
elif [ $error -eq 1 ]
then
   exit 1
else
   exit 1
fi

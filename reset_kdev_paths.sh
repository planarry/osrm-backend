prname=`pwd`"/.kdev4/"`ls ./.kdev4`
dst=`pwd`
src=$( grep -oP '(?<=Build Directory Path=file://).*(?=/build)' "$prname" )
echo replace $src on $dst
src=$( echo $src | sed 's/\//\\\//g' )
dst=$( echo $dst | sed 's/\//\\\//g' )
sed -i 's/'$src'/'$dst'/g' "$prname"
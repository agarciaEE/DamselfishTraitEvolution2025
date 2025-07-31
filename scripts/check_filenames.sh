while read -r name; do
# Ignore comment lines
   if [[ -f "$2/$name.jpg" ]]
   then
      echo  "$name.jpg exists"
   else
      echo "ERROR!! No File for $name.jpg"
      exit 1
   fi
done < "$1" 
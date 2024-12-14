BEGIN        {  printf("i pkginfo=./pkginfo\n") }
$1 ~ /^[fd]/ {  if ($3 != "out")
	          print $1, $2, $3, $4, "bin", "bin"
             } 

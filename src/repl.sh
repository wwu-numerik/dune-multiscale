while read data; do
	result=$data

#some color declarations
	GREEN=`echo -e '\033[49;32m'`
	NORMAL=`echo -e '\033[0m'`
	RED=`echo -e '\033[49;31m'`
	BROWN=`echo -e '\033[49;33m'`
	BLUE=`echo -e '\033[49;34m'`
	PINK=`echo -e '\033[49;35m'`
	CYAN=`echo -e '\033[49;36m'`
	WHITE=`echo -e '\033[49;1;37m'`
	LIGHTGRAY=`echo -e '\033[49;0;37m'`
	ULINE=`echo -e '\033[49;4;30m'`
	REDINVERTED=`echo -e '\033[49;7;31m'`
	BROWNINVERTED=`echo -e '\033[49;7;33m'`
	BLACK=`echo -e '\033[49;30m'`
	BLACK=$WHITE

#some colorings
#change it - if needed
#change with care
	COLORTEMPLATEBRACKET1=`echo -e '\033[38;5;21m'`
	COLORTEMPLATEBRACKET2=`echo -e '\033[38;5;93m'`
	COLORTEMPLATEBRACKET3=`echo -e '\033[38;5;129m'`
	COLORTEMPLATEBRACKET4=`echo -e '\033[38;5;165m'`
	COLORTEMPLATEBRACKET5=`echo -e '\033[38;5;162m'`
	COLORTEMPLATEBRACKET6=`echo -e '\033[38;5;160m'`
	COLORTEMPLATEBRACKET7=`echo -e '\033[38;5;166m'`
	COLORTEMPLATEBRACKET8=`echo -e '\033[38;5;172m'`
	COLORTEMPLATEBRACKET9=`echo -e '\033[38;5;176m'`
	COLORTEMPLATEBRACKET10=`echo -e '\033[38;5;184m'`
	COLORTEMPLATEBRACKET11=`echo -e '\033[38;5;1m'`
	COLORTEMPLATEBRACKET12=`echo -e '\033[38;5;3m'`
	COLORTEMPLATEBRACKET13=`echo -e '\033[38;5;29m'`
	COLORTEMPLATEBRACKET14=`echo -e '\033[38;5;65m'`
	COLORTEMPLATEBRACKET15=`echo -e '\033[38;5;62m'`
	COLORTEMPLATEBRACKET16=`echo -e '\033[38;5;60m'`
	COLORTEMPLATEBRACKET17=`echo -e '\033[38;5;66m'`
	COLORTEMPLATEBRACKET18=`echo -e '\033[38;5;72m'`
	COLORTEMPLATEBRACKET19=`echo -e '\033[38;5;78m'`
	COLORTEMPLATEBRACKET20=`echo -e '\033[38;5;84m'`

	COLORFILEPATH=$LIGHTGRAY
	COLORTEMPLATEREPLACE=$BLUE	
	COLORERROR=$REDINVERTED
	COLORWARNING=$BROWNINVERTED
	COLORIMPORTANTWORD=$BLUE

	REMOVEPATH="\/u\/s_girk01\/"



echo -e "$result" | sed -e "s/\([^\-]\)\([<>]\)/\1 \2 /g" | awk '{
												auf = 0
												Anz = split($0, array, " ")
												for(i=1; i <= Anz; i++){
													if ( match (array[i],"<") ){
														printf "%s%i ", array[i], auf
														auf++;
													}
													else if( match (array[i],">") &&  !match (array[i],"->") ){
														auf--
														printf "%s%i ", array[i], auf																
													}
													else
														printf "%s ", array[i]														
												}
												printf "\n"
												}' | sed -e  "s/<0\( [^<]*\)<1/$COLORTEMPLATEBRACKET1<\1$BLACK<1/g;\
												s/>1\( [^<]*\)<1/>1$COLORTEMPLATEBRACKET1\1$BLACK<1/g;\
												s/>1\( [^<]*\)>0/>1$COLORTEMPLATEBRACKET1\1>$BLACK/g;\
												s/<0\( [^<]*\)>0/$COLORTEMPLATEBRACKET1<\1>$BLACK/g;\
												s/<1\( [^<]*\)<2/$COLORTEMPLATEBRACKET2<\1$BLACK<2/g;\
												s/>2\( [^<]*\)<2/>2$COLORTEMPLATEBRACKET2\1$BLACK<2/g;\
												s/>2\( [^<]*\)>1/>2$COLORTEMPLATEBRACKET2\1>$BLACK/g;\
												s/<1\( [^<]*\)>1/$COLORTEMPLATEBRACKET2<\1>$BLACK/g;\
												s/<2\( [^<]*\)<3/$COLORTEMPLATEBRACKET3<\1$BLACK<3/g;\
												s/>3\( [^<]*\)<3/>3$COLORTEMPLATEBRACKET3\1$BLACK<3/g;\
												s/>3\( [^<]*\)>2/>3$COLORTEMPLATEBRACKET3\1>$BLACK/g;\
												s/<2\( [^<]*\)>2/$COLORTEMPLATEBRACKET3<\1>$BLACK/g;\
												s/<3\( [^<]*\)<4/$COLORTEMPLATEBRACKET4<\1$BLACK<4/g;\
												s/>4\( [^<]*\)<4/>4$COLORTEMPLATEBRACKET4\1$BLACK<4/g;\
												s/>4\( [^<]*\)>3/>4$COLORTEMPLATEBRACKET4\1>$BLACK/g;\
												s/<3\( [^<]*\)>3/$COLORTEMPLATEBRACKET4<\1>$BLACK/g;\
												s/<4\( [^<]*\)<5/$COLORTEMPLATEBRACKET5<\1$BLACK<5/g;\
												s/>5\( [^<]*\)<5/>5$COLORTEMPLATEBRACKET5\1$BLACK<5/g;\
												s/>5\( [^<]*\)>4/>5$COLORTEMPLATEBRACKET5\1>$BLACK/g;\
												s/<4\( [^<]*\)>4/$COLORTEMPLATEBRACKET5<\1>$BLACK/g;\
												s/<5\( [^<]*\)<6/$COLORTEMPLATEBRACKET6<\1$BLACK<6/g;\
												s/>6\( [^<]*\)<6/>6$COLORTEMPLATEBRACKET6\1$BLACK<6/g;\
												s/>6\( [^<]*\)>5/>6$COLORTEMPLATEBRACKET6\1>$BLACK/g;\
												s/<5\( [^<]*\)>5/$COLORTEMPLATEBRACKET6<\1>$BLACK/g;\
												s/<6\( [^<]*\)<7/$COLORTEMPLATEBRACKET7<\1$BLACK<7/g;\
												s/>7\( [^<]*\)<7/>7$COLORTEMPLATEBRACKET7\1$BLACK<7/g;\
												s/>7\( [^<]*\)>6/>7$COLORTEMPLATEBRACKET7\1>$BLACK/g;\
												s/<6\( [^<]*\)>6/$COLORTEMPLATEBRACKET7<\1>$BLACK/g;\
												s/<7\( [^<]*\)<8/$COLORTEMPLATEBRACKET8<\1$BLACK<8/g;\
												s/>8\( [^<]*\)<8/>8$COLORTEMPLATEBRACKET8\1$BLACK<8/g;\
												s/>8\( [^<]*\)>7/>8$COLORTEMPLATEBRACKET8\1>$BLACK/g;\
												s/<7\( [^<]*\)>7/$COLORTEMPLATEBRACKET8<\1>$BLACK/g;\
												s/<8\( [^<]*\)<9/$COLORTEMPLATEBRACKET9<\1$BLACK<9/g;\
												s/>9\( [^<]*\)<9/>9$COLORTEMPLATEBRACKET9\1$BLACK<9/g;\
												s/>9\( [^<]*\)>8/>9$COLORTEMPLATEBRACKET9\1>$BLACK/g;\
												s/<8\( [^<]*\)>8/$COLORTEMPLATEBRACKET9<\1>$BLACK/g;\
                        s/<9\( [^<]*\)<10/$COLORTEMPLATEBRACKET10<\1$BLACK<10/g;\
												s/>10\([^<]*\)<10/>10$COLORTEMPLATEBRACKET10\1$BLACK<10/g;\
												s/>10\([^<]*\)>9/>10$COLORTEMPLATEBRACKET10\1>$BLACK/g;\
												s/<9\( [^<]*\)>9/$COLORTEMPLATEBRACKET10<\1>$BLACK/g;\
                        s/<10\([^<]*\)<11/$COLORTEMPLATEBRACKET11<\1$BLACK<11/g;\
												s/>11\([^<]*\)<11/>11$COLORTEMPLATEBRACKET11\1$BLACK<11/g;\
												s/>11\([^<]*\)>10/>11$COLORTEMPLATEBRACKET11\1>$BLACK/g;\
												s/<10\([^<]*\)>10/$COLORTEMPLATEBRACKET11<\1>$BLACK/g;\
												s/<11\([^<]*\)<12/$COLORTEMPLATEBRACKET12<\1$BLACK<12/g;\
												s/>12\([^<]*\)<12/>12$COLORTEMPLATEBRACKET12\1$BLACK<12/g;\
												s/>12\([^<]*\)>11/>12$COLORTEMPLATEBRACKET12\1>$BLACK/g;\
												s/<11\([^<]*\)>11/$COLORTEMPLATEBRACKET12<\1>$BLACK/g;\
												s/<12\([^<]*\)<13/$COLORTEMPLATEBRACKET13<\1$BLACK<13/g;\
												s/>13\([^<]*\)<13/>13$COLORTEMPLATEBRACKET13\1$BLACK<13/g;\
												s/>13\([^<]*\)>12/>13$COLORTEMPLATEBRACKET13\1>$BLACK/g;\
												s/<12\([^<]*\)>12/$COLORTEMPLATEBRACKET13<\1>$BLACK/g;\
												s/<13\([^<]*\)<14/$COLORTEMPLATEBRACKET14<\1$BLACK<14/g;\
												s/>14\([^<]*\)<14/>14$COLORTEMPLATEBRACKET14\1$BLACK<14/g;\
												s/>14\([^<]*\)>13/>14$COLORTEMPLATEBRACKET14\1>$BLACK/g;\
												s/<13\([^<]*\)>13/$COLORTEMPLATEBRACKET14<\1>$BLACK/g;\
												s/<14\([^<]*\)<15/$COLORTEMPLATEBRACKET15<\1$BLACK<15/g;\
												s/>15\([^<]*\)<15/>15$COLORTEMPLATEBRACKET15\1$BLACK<15/g;\
												s/>15\([^<]*\)>14/>15$COLORTEMPLATEBRACKET15\1>$BLACK/g;\
												s/<14\([^<]*\)>14/$COLORTEMPLATEBRACKET15<\1>$BLACK/g;\
												s/<15\([^<]*\)<16/$COLORTEMPLATEBRACKET16<\1$BLACK<16/g;\
												s/>16\([^<]*\)<16/>16$COLORTEMPLATEBRACKET16\1$BLACK<16/g;\
												s/>16\([^<]*\)>15/>16$COLORTEMPLATEBRACKET16\1>$BLACK/g;\
												s/<15\([^<]*\)>15/$COLORTEMPLATEBRACKET16<\1>$BLACK/g;\
												s/<16\([^<]*\)<17/$COLORTEMPLATEBRACKET17<\1$BLACK<17/g;\
												s/>17\([^<]*\)<17/>17$COLORTEMPLATEBRACKET17\1$BLACK<17/g;\
												s/>17\([^<]*\)>16/>17$COLORTEMPLATEBRACKET17\1>$BLACK/g;\
												s/<16\([^<]*\)>16/$COLORTEMPLATEBRACKET17<\1>$BLACK/g;\
												s/<17\([^<]*\)<18/$COLORTEMPLATEBRACKET18<\1$BLACK<18/g;\
												s/>18\([^<]*\)<18/>18$COLORTEMPLATEBRACKET18\1$BLACK<18/g;\
												s/>18\([^<]*\)>17/>18$COLORTEMPLATEBRACKET18\1>$BLACK/g;\
												s/<17\([^<]*\)>17/$COLORTEMPLATEBRACKET18<\1>$BLACK/g;\
												s/<18\([^<]*\)<19/$COLORTEMPLATEBRACKET19<\1$BLACK<19/g;\
												s/>19\([^<]*\)<19/>19$COLORTEMPLATEBRACKET19\1$BLACK<19/g;\
												s/>19\([^<]*\)>18/>19$COLORTEMPLATEBRACKET19\1>$BLACK/g;\
												s/<18\([^<]*\)>18/$COLORTEMPLATEBRACKET19<\1>$BLACK/g;\
												s/<-[0-9]*/</g;\
												s/>-[0-9]*/>/g;\
												s/<[0-9]*/</g;\
												s/>[0-9]*/>/g;\
												s/$REMOVEPATH//;\
												s/\[with /\n\[with /g;\
												s/\(^.[^:]*:[0-9]*: \)/$COLORFILEPATH\1$BLACK\n\t/i;\
												s/\(^.[^:]*:[0-9]*:[0-9]*: \)/$COLORFILEPATH\1$BLACK\n\t/i;\
												s/\(^.[^:]*:\)\([^0-9]\)/$COLORFILEPATH\1$BLACK\n\t\2/;\
												s/\(-[DIM]\)/$COLORIMPORTANTWORD\1$BLACK/g;\
												s/»\([^«]*\)«/\n\t$ULINE\1$NORMAL/g;\
											 	s/= \([A-Za-z0-9_:.~]*\)/= $COLORTEMPLATEREPLACE\1$BLACK/g;\
												s/\(error\)\([^a-z]\)/$COLORERROR\1$NORMAL\2/ig;\
												s/\(fehler\)\([^a-z]\)/$COLORERROR\1$NORMAL\2/ig;\
												s/\(warn[iu]ng\)\([^a-z]\)/$COLORWARNING\1$NORMAL\2/ig;\
												s/make/$COLORIMPORTANTWORD&$BLACK/g;\
												s/configure/$COLORIMPORTANTWORD&$BLACK/g;\
												s/*/\*/g"

											 	#s/= \([^, ]*\),/= $COLORTEMPLATEREPLACE\1$BLACK,/g;\
done


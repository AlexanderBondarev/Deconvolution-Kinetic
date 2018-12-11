#include <stdio.h>
#include <stdlib.h>

int main()
{
 int flag=0;
 char s[4096];

 while(!feof(stdin))
  {
    gets(s);
    if(strstr(s,"Signal correct")!=NULL) flag=1;
    if(strstr(s,"Time")!=NULL) printf("%s\n",s);
    if(flag==1) printf("%s\n",s);
    if(strstr(s,"Ampoule removed")!=NULL) break;
  }
}

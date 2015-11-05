int n2Z(char *elem){
  if(strcmp(elem, "H")==0) return 1;
  else if (strcmp(elem, "He")==0) return 2;
  else if (strcmp(elem, "Li")==0) return 3;
  else if (strcmp(elem, "Be")==0) return 4;
  else if (strcmp(elem, "B" )==0) return 5;
  else if (strcmp(elem, "C" )==0) return 6;
  else if (strcmp(elem, "N" )==0) return 7;
  else if (strcmp(elem, "O" )==0) return 8;
  else if (strcmp(elem, "F" )==0) return 9;
  else if (strcmp(elem, "Ne")==0) return 10;
  else if (strcmp(elem, "Na")==0) return 11;
  else if (strcmp(elem, "Mg")==0) return 12;
  else if (strcmp(elem, "Al")==0) return 13;
  else if (strcmp(elem, "Si")==0) return 14;
  else if (strcmp(elem, "P" )==0) return 15;
  else if (strcmp(elem, "S" )==0) return 16;
  else if (strcmp(elem, "Cl")==0) return 17;
  else if (strcmp(elem, "Ar")==0) return 18;
  else if (strcmp(elem, "K" )==0) return 19;
  else if (strcmp(elem, "Ca")==0) return 20;
  else{
    printf("ERROR! element %s not found!\n", elem);
    exit(1);
  }
}

double sqr(double x){return x*x;}

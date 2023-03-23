#include <stdio.h> 
#include <math.h> /* need this for powers */
#include <stdlib.h>

// program to calculate binomial probabilities

double factorial(double number) 
{
  if(number>1) return number * factorial(number-1);
  else return 1;
} 

// THESE DOUBLE WERE ORIGINAL INTERVALS, BUT ONCE NUMBERS GET LARGE GOES PEAR SHAPED
// SEE http://www.aspire.cs.uah.edu/textbook/CPP7002.html
main()
{
  float  p, q;
  long double  total_p, total_q, k_fac, n_fac,d_fac,m_fac,K_fac,kay_fac, diff_fac, E, V; 
  long double *P_k, pk,total_pk,Z, Z_b;
  int i, K, n, k, d, number, m, diff;
  int index = 0, j = 0;
  long double *P, *Q, *e, *v; //, *mean;
  char systemCall[100];

  printf("\nProgram determines binary [binomial] probability\nEnter the probability of success (0.5 for even odds):"); 
  scanf("%f", &p);

  while (p >1){
    printf("Illegal value, try again!\n");
    scanf("%f", &p);
  }

  q = 1-p;  //probability of failure

  printf("\nNumber of trials [e.g. coin tosses]: "); 
  scanf("%d", &n); 

  printf("Number of positive outcomes [e.g. heads]: "); 
  scanf("%d", &k);

  P = (long double*)malloc(n*sizeof(long double));  //allocating memory which is at least n
  e = (long double*)malloc(n*sizeof(long double));   
  v = (long double*)malloc(n*sizeof(long double));
 
  while (k>n){
    printf("This cannot be more than the number of trials knucklehead! Enter again: \n");
    scanf("%d", &k);
  } 

  //printf("k = %d n = %d p = %1.1e np = %1.2e\n", k, n, p, (float)n*p);  
  Z = ((float)k+0.5 - p*n)/pow(p*n*(1-p),0.5);
  //  Z_b = ((float)k-0.5 - p*n)/pow(p*n*(1-p),0.5);//http://courses.wcupa.edu/rbove/Berenson/10th%20ed%20CD-ROM%20topics/section6_5.pdf

  K = n -k;  //printf("n = %d k = %d, K = %d\n",n, k, K);
  
  for (k=k; k<=n; ++k){ 
    for (number =1; number<=k; ++number); //printf("%d!= %f\n", number, factorial(number)); 
    k_fac = factorial(number-1);  // printf("k_fac =%f\n",k_fac);
     
    for (number =1; number <= n; ++number);// printf("%d!= %f\n", number, factorial(number)); 
    n_fac = factorial(number-1); // printf("n_fac =%f\n",n_fac);
  
    //////////////////////////////////////////////////////////////////////
    d=n-k; 
    if (d >0){
      for (number =1; number<=d; ++number); d_fac = factorial(number-1); 
    }
    else d_fac = 1;
    /////////////////////////////////////////////////////////////////////
    //P(k out of n trials) = p^{k} * q^{n-k} * n! / k!(n-k)!
    P[index] = (pow(p,k))*(pow(q,n-k))*n_fac/(k_fac*d_fac);
    // printf("Probability of %d 'heads' out of %d 'tosses' is %e\n", k, n, P[index]);
    index++;     
  }
  total_p = 0; 
  for(i=0; i<n; i++) { 
    total_p += P[i];
  }

  printf("-------------------------------------------------------------\n");
  printf("For p = %1.3f and q = %1.3f ...\n", p, q); 
  printf("Probability of >= %d [or more] out of %d is then %1.4Le\n", k - index, n, total_p);
   
  int kay = k-index; // to save this number as gets loast

  ///////////////////////////// POISSON APPROX //////////////////////////
  P_k = (long double*)malloc(n*sizeof(long double));
 
  for (number =1; number<=kay; ++number)// printf("%d!= %f\n", number, factorial(number)); 
    kay_fac = factorial(number);
  // printf("n = %d k = %d, K = %d kay = %d kay! = %1.2f\n",n, k, K, kay, kay_fac);
  int pois = kay;
  P_k[index] = exp(-1*n*p)*pow(n*p,(kay))/kay_fac;

  ////////////////////////// NEW BIT TO GET FEWER THAN ALSO /////////////
  i = index = 0; 
  Q = (long double*)malloc(n*sizeof(long double)); // or Bus error
    
  for (K=K; K<=n; ++K){ 
    for (number =1; number<=K; ++number); //printf("%d!= %f\n", number, factorial(number)); 
    K_fac = factorial(number-1);  // printf("K_fac =%f\n",K_fac);
     
    for (number =1; number <= n; ++number);// printf("%d!= %f\n", number, factorial(number)); 
    n_fac = factorial(number-1); // printf("n_fac =%f\n",n_fac);
  
    d=n-K; 
    if (d >0){
      for (number =1; number<=d; ++number); d_fac = factorial(number-1); 
    }
    else d_fac = 1;
     
    // printf("K = %d kay = %d\n", K, kay);
    if (kay > 0) Q[index] = (pow(q,K))*(pow(p,n-K))*n_fac/(K_fac*d_fac);  // note p and q swapped around - using kay at lost K
    else {
      Q[index] = pow(q,n); //Q[index] = (pow(p,K))*(pow(q,n-K)); // K = 0 so K! =1 so K!(n-K)! =  n!  so n_facs cancel 
      printf("Probability of 0 out of %d is then %1.4Le\n",  n, Q[index]);
    }
    index++;     
  }
   
  total_q = 0; 
  if (kay > 0) for(i=0; i<n; i++) total_q += Q[i];
  else  total_q = Q[index];

  if (kay > 0) printf("Probability of <= %d [or fewer] out of %d is then %1.4Le\n", kay, n, total_q);

  printf("\n---------------- NORMAL APPROXIMATION ---------------------\n");

  /* A continuity correction factor is used when you use a continuous function to approximate a discrete one. For example, when you want to approximate a binomial with a normal distribution. According to the Central Limit Theorem, the sample mean of a distribution becomes approximately normal if the sample size is large enough. The binomial distribution can be approximated with a normal distribution too, as long as n*p and n*q are both at least 5. */

  /* The continuity correction factor a way to account for the fact that a normal distribution is continuous, and a binomial is not. When you use a normal distribution to approximate a binomial distribution, you’re going to have to use a continuity correction factor. It’s as simple as adding or subtracting .5 to the discrete x-value: use the following table to decide whether to add or subtract. */

  //http://www.statisticshowto.com/what-is-the-continuity-correction-factor/
 
  float mean, SD, zvalue_leq, zvalue_geq, zvalue_gt, zvalue_lt;

  mean = n*p; SD = sqrt(mean*q); zvalue_leq = zvalue_gt = (kay+0.5 - mean)/SD; zvalue_geq = zvalue_lt = (kay-0.5 - mean)/SD;

  //printf("k = %d\n", kay); printf("   Mean is %1.5g and SD = %1.5g\n", mean, SD);
   
  ///////  INV SIGMA /////////////////////////
  int npts = 99999;
  double prob, whole_prob,x[npts], y[npts], A, a[npts],*area; //, Y[npts], *y;
  float zvalue, top =20; // up to 20 simga
  float norm = pow((2*3.141592654),0.5); // normalisation

  if (zvalue_geq <0)  zvalue = -1*zvalue_geq;
  else zvalue = zvalue_geq;

  for (i=0; i<npts; ++i){ 
    a[0] = 0; // first point scres everything up
    x[i] = zvalue + ((top-zvalue)/npts)*i; 
    y[i]=(1/norm)*exp(-0.5*x[i]*x[i]); 

    a[i] = (x[i]-x[i-1])*y[i];  
    a[i]=a[i]+a[i-1];  // summing up
    area = &a[i];
    A = *area;
  }
  //  printf("Probability of >= %d [or more] out of %d gives Z-value = %1.3f => P = %1.3e\n", kay, n, zvalue_geq, 1-A);
  //printf("      of < %d out of %d gives Z-value = %1.3f => P = %1.3e\n", kay, n, zvalue_geq, A);

  //////////////////////////////////////////////////////// 
  if (zvalue_leq <0)  zvalue = -1*zvalue_leq;
  else zvalue = zvalue_leq;

  for (i=0; i<npts; ++i){ 
    a[0] = 0; // first point scres everything up
    x[i] = zvalue + ((top-zvalue)/npts)*i; 
    y[i]=(1/norm)*exp(-0.5*x[i]*x[i]); 

    a[i] = (x[i]-x[i-1])*y[i];  
    a[i]=a[i]+a[i-1];  // summing up
    area = &a[i];
    A = *area;
  }
  printf("Probability of >= %d [or more] out of %d gives Z-value = %1.3f => P = %1.3e\n", kay, n, zvalue, 1-A);
  printf("Probability of <= %d [or fewer]  out of %d gives Z-value = %1.3f => P = %1.3e\n", kay, n, zvalue, A);
  
  
 
  // printf("\n For scripts...\n%1.3f \n%1.3e\n", zvalue, A);

  ///////////////////////////////////////////////////////////////////////////////////////////////////////

  /*  printf("          Using normal approx. Z value of %d or fewer of %d \n             is %1.3Lf sigma [run inv_sigma on this?]\n", kay, n, Z); //as k now used */
  /*     ///////////////////////////////// MEAN ////////////////////////////////////      */
  /*     for (m=1; m<=n; ++m){ */
  /*       for (number =1; number<=m; ++number){ */
  /* 	m_fac = factorial(number); */
  /*       } */
  /*       for (number =1; number <= n; ++number){ */
  /* 	n_fac = factorial(number); */
  /*       } */
  /*       d=n-m;  //redefine d  */
  /*       if (d >0){ */
  /* 	for (number =1; number<=d; ++number);{ */
  /* 	  d_fac = factorial(number-1);  //don't know why this is different - just works */
  /* 	  // printf("d_fac =%f\n",d_fac); */
  /* 	} */
  /*       } */
  /*       else d_fac = 1; */
     
  /*       e[index] =  m * (pow(p,m))*(pow(q,n-m))*n_fac/(m_fac*d_fac); */
  /*       index++; */
  /*       //have to tally this up and close loop */
  /*       e[index+1]= e[index] + e[index-1];  */
  /*       //printf("d_fac =%1.1f m = %d, E = %1.1f e = %1.4e\n",d_fac, m, E,e[index]); */
  /*     } */
   
  /*     mean = &e[index+1];    // BELONGS TO # */
  /*     E = *mean;          // # */
     
  /*     //  printf("\nMean = %1.3f [cf. np = %1.3f, as used by http://faculty.vassar.edu/lowry/binomialX.html]", E, n*p);  */
  
  /*     v = (long double*)malloc(n*sizeof(long double)); // HAD TO RESTATE OR FOLLOWING RUNS OUT OF STEAM!! */
  /*     ///////////////////////VARAINCE AND SD /////////////////////////////      */
  /*     for (m=1; m<=n; ++m){ */
  /*       for (number =1; number<=m; ++number){ */
  /* 	m_fac = factorial(number); */
  /*       } */
  /*       for (number =1; number <= n; ++number){ */
  /* 	n_fac = factorial(number); */
  /*       } */
  /*       d=n-m; // redefine varinnce */
  /*       if (d >0){ */
  /* 	for (number =1; number<=d; ++number);{ */
  /* 	  d_fac = factorial(number -1);  //don't know why this is different - just works */
  /* 	  // printf("number = %d, m = %d, d_fac =%1.0g\n",number, m, d_fac); */
  /* 	} */
  /*       } */
  /*       else  d_fac = 1; */
          
  /*       v[index] =  pow((float)m-E,2) * (pow(p,m))*(pow(q,n-m))*n_fac/(m_fac*d_fac); */
  /*        index++; */
  /*        v[index+1]= v[index] + v[index-1];  */
  /*        // printf("number = %d, d_fac =%1.1e m = %d, n = %d, E = %1.1f var = %1.4e\n",number, d_fac, m, n, E,v[index]); */
  /*     } */
   
  /*     // printf("Variance = %1.3f\n", v[index+1]);  */
  /*     printf("\nMean = %1.3e +/- %1.3e [sigma]\n", n*p, sqrt(v[index+1]));  */
  /*     printf("-----------------------------------------------------------------\n"); */
  
  // PRINT INPUT VALUES TO SEE WHAT MONTE CARLO IN PYTHON GIVES
  // printf("n = %d k = %d", n, k);
  //sprintf(systemCall,"/Users/stephencurran/python/./bin-prob.py %d %d %f",n, kay, p); system(systemCall); 
    
}

//                cc -o bin-prob bin-prob.c -lm ; ./bin-prob




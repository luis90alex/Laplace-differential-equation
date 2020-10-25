
#include "Myheader.h"

void previo(double *L,double *dx,double *a, int *N, int *k2,int *k1){

    *N=(int) (*L/(*dx));   /// el nombre de punts és una variable entera
    *k2=(int) (*a/(*dx));

    /// per fer la malla simètrica s'han de cambiar lleugerament els paràmetres inicials
    ///  N=k1+k2+k1 ; per a que hi hagui el mateix nombre de punts entre 0 fins k1  i  des de (k1+k2) fins a N, s'han d'imposar les condicions de sota.

    if ( ( *N- *k2)%2 == 0 ) {
		*k1 = ( *N- *k2)/2;
	}
	else {
		*N = *N-1;
		*k1 = (*N-*k2)/2;
	}
  printf("N es %i\n",*N);/// s'imprimeix per pantalla el valor final de N

}


void introduce(int N,int n,int k1,int k2,int A,double **m1){
    int i,j,k,l;/// els contadors que s'utilitzaràn


    for(i=0;i<=N;i++){  ///primer s'ompliran les primeres i les últimes N files de la matriu m1
            l=i+N*(N+1);/// aquestes son les últimes N files
		for(j=0;j<n+1;j++){/// cada j representa la columna de m1, j arriba fins a n atès que s'ha inclòs el vector b en la matriu m1

			if(j==i){/// s'imposa l'equació de laplace per als punts de la vora de sota de la malla
                m1[i][j]=4;
            }
            else {
                 if(j==i+1 || j==i-1 || j==i+N+1 ){
                       m1[i][j]=-1;
                 }
                 else {
                       m1[i][j]=0;
                 }
            }

            if(j==l){/// s'imposa l'equació de laplace per als punts de la vora de dalt de la malla
                m1[l][j]=4;
            }
            else {
                 if(j==l+1 || j==l-1 || j==l-(N+1) ){
                       m1[l][j]=-1;
                 }
                 else {
                       m1[l][j]=0;
                 }
            }
		}
	}
	for(k=1;k<=N-1;k++){ ///aquí s'omplen les files de la matriu que corresponen als punts de la vora esquerra i dreta de la malla
            i=k*(N+1);l=N+k*(N+1);
            /// amb i es designaran les files de la vora esquerra i amb l les de la vora dreta,s'han omès les files que ja han sigut utilitzades
            /// és per això que k comença en 1 i acaba en N-1
        for(j=0;j<n+1;j++){

			if(j==i){ /// s'imposa l'equació de laplace per als punts de la vora esquerra de la malla
                m1[i][j]=4;
            }
            else {
                 if(j==i+1 || j==i-(N+1) || j==i+N+1 ){
                       m1[i][j]=-1;
                 }
                 else {
                       m1[i][j]=0;
                 }
            }

            if(j==l){ /// s'imposa l'equació de laplace per als punts de la vora dreta de la malla
                m1[l][j]=4;
            }
            else {
                 if(j==l-1 || j==l-(N+1) || j==l+(N+1) ){
                       m1[l][j]=-1;
                 }
                 else {
                       m1[l][j]=0;
                 }
            }

        }
	}

	///Aquí corregim tres números que per com hem fet el bucle posaven de manera incorrecta i que no tocaven

	//cantonada inferior dreta
	m1[N][N+1]=0;

	//cantonada superior esquerra
	m1[N*(N+1)][N*(N+1)-1] = 0;

	//cantonada superior dreta
	m1[N*(N+1)+N][N*(N+1)+N+1] = 0;


	///aquí s'omplen les files corresponents als punts que no es situen a cap vora de la malla
	for(k=1;k<=N-1;k++){
        for(l=1;l<=N-1;l++){
            i=l+k*(N+1);/// es designa la fila i-èsima de la matriu mitjançant l i k, no estan incloses les files de les vores
            for(j=0;j<=n;j++){
                if( k>=k1-1 && k<=k1+k2+1 && l>=k1-1 && l<=k1+k2+1 ){/// aquesta condició és per a tots els punts del quadrat interior i del quadrat "ampliat"

                    if(j==n){/// s'imposa que els punts valguin A
                         m1[i][j]=A;
                    }
                    else {
                         m1[i][j]=(i==j)? 1:0;
                    }
                }
                else {///condició per els punts que no son a les vores ni al quadrat "ampliat"
                        if(j==i){
                            m1[i][j]=4;
                        }
                        else {
                            if(j==i+1 || j==i-1 || j==i+N+1 || j==i-N-1){
                                m1[i][j]=-1;
                            }
                            else {
                                m1[i][j]=0;
                            }
                        }

                } ///aquí es on fem la petita correcció  a la condició de contorn de les cares

                        if( (k==k1-1 && l==k1-1) || (k==k1-1 && l==k1+k2+1) || (k==k1+k2+1 && l==k1-1) || (k==k1+k2+1 && l==k1+k2+1) ){
                            if(j==i){
                                m1[i][j]=4;
                            }
                            else {
                                if(j==i+1 || j==i-1 || j==i+N+1 || j==i-N-1){
                                    m1[i][j]=-1;
                                }
                                else {
                                    m1[i][j]=0;
                                }
                            }

                        }
            }
        }
	}

}

void triangular(int N,int n,double **m1){///funció que convertirà m1 en triangular superior

    int i,j,k;
	double M;

	for(k=0;k<n;k++){/// a mesura que k augmenta, la columna (k-1)-èsima de la matriu té zeros per sota de la diagonal
		for(i=k+1;i<n;i++){
            M=(m1[i][k])/(m1[k][k]);/// aquest és el coeficient que s'utilitzarà per fer zeros per sota de la diagonal amb gauss
			for(j=0;j<n+1;j++){

			if(i<=k){
				m1[i][j]=m1[i][j];
			}
			else{
				m1[i][j]=(m1[i][j]-M*m1[k][j]);
			}
			}
		}

	}
}

void solution(int N,int n,double dx, double **m1, double *sol){///aqui es solucionara el sistema per substitució regresiva

    int i,j,s;
	double sum=0;
    sol[n-1]=m1[n-1][n]/(m1[n-1][n-1]); ///és el valor de l'últim punt de la malla, el de la cantonada superior dreta
    /// es comença amb aquesta solució i s'obtenen les altres a partir de l'última

	for(i=n-2;i>=0;i--){
		for(j=i+1;j<=n-1;j++){
			sum=sum+m1[i][j]*sol[j];
		}
	sol[i]=(m1[i][n]-sum)/m1[i][i];
	sum=0;
	}

	FILE *fichero;/// es disposen totes les solucions en un arxiu .txt per tal de poder-lo representar posteriorment
	char laplace[]="solucioneslaplace.txt";/// el nom de l'arxiu
	fichero=fopen(laplace,"w");/// s'obre l'arxiu per escriure-hi en ell

	for(i=0;i<=N;i++){
            for(j=0;j<=N;j++){
                s=(N+1)*j+i;/// la variable s relaciona les solucions de la funció a cada punt amb aquell punt
                fprintf(fichero,"%lf %lf %lf \n",i*dx,j*dx,sol[s]);
            }
	}

}


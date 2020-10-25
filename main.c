
#include "Myheader.h"

int main()
{
    double L,dx,a,A,**m1,*sol;/// Es declaren 6 variables on L �s la longitud de la malla,dx ser� l'increment tant en x com en y
    /// a �s la longitud del cuadrat interior on el potencial val A, m1 ser� la matriu que es forma per resoldre el sistema
    /// sol ser� el vector que cadascuna de les seves components  inclour� el valor de la funci� en un punt determinat

    int N,n,k1,k2;/// N representa el nombre de punts de la malla,n �s el nombre de files  de la matriu m1 (n x(n+1);
    /// k1 �s el nombre de punts que van desde 0 fins quan comen�a el cuadrat interior tal com s'ha explicat a l'informe
    /// k2 �s el nombre de punts que t� una de les cares del cuadrat interior

    char myname[]="Luis Farfan Moina i Esteban Jeronimo Piedrahita";/// varible on emmagatzemmar els nostres noms

    printf("Fet per  %s\n",myname);
    printf("Introdueix L,dx,a i A (en aquest ordre)\n");/// es demana per pantalla que s'introdueixon els valors de L,dx,a,A
    scanf("%lf %lf %lf %lf",&L,&dx,&a,&A);/// es recullen aquests valors en les variables definides anteriorment


    previo(&L,&dx,&a,&N,&k2,&k1);/// funci� que cambiar� lleugerament els par�metres inicials per tal de que el programa funcioni correctament
    n= (N)*(N)+2*(N)+1;


    m1 = (double**)calloc(n,sizeof(double));/// declaracio de la matriu m1
    int i;
    for (i=0 ; i <= n ; i++) {///bucle per inicialitzar la matriu m1
		m1 [i] = (double*)calloc(n,sizeof(double));
	}
	sol = (double*)calloc(n,sizeof(double));///declaraci� del vector sol


    introduce(N,n,k1,k2,A,m1);/// funci� que creara la matriu m1
    triangular(N,n,m1);/// funci� que transformara la matriu m1 en una matriu triangular mitjan�ant gauss
    solution(N,n,dx,m1,sol);/// funci� que resol la matriu m1 triangular



   return 0;
}

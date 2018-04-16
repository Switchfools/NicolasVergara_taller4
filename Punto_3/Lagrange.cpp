//
//  Lagrange.cpp
//  
//
//  Created by Nicolas Felipe Vergara Duran on 15/04/18.
//
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <complex>
#include <math.h>
using namespace std;
int Number_data(){
    string line;
    ifstream myfile ("datos.txt");
    if (myfile.is_open())
    {
        int numero_datos=0;
        while ( getline (myfile,line) )
        {
            numero_datos++;
        }
        return (numero_datos);
        myfile.close();
    }
    
    else cout << "Unable to open file";
}
double** DatatoMatrix(int Filas, int Columnas){
    double **Matrix = new double*[Filas];
    for(int i=0;i<Filas;i++){
        Matrix[i] = new double[Columnas];
    }
    string line;
    ifstream myfile ("First.txt");
    if (myfile.is_open())
    {
        int i=0;
        while ( getline (myfile,line) )
        {
             myfile >> Matrix[i][0] >> Matrix[i][1];
            i++;
        }
        myfile.close();
    }
    return Matrix;
}
double* linspace(double in, double end, int size){
    double* vector = new double[size];
    double h=(end-in)/size;
    double point=in;
    for(int i=0;i<size;i++){
        vector[i]=point;
        point+=h;
    }
    return vector;
}
double** PolinomiodeLagrange(int n_datos,double** datos){
    double* new_x=linspace(datos[0][0],datos[n_datos-2][0],n_datos);
    double **Matrix = new double*[n_datos];
    for(int i=0;i<n_datos;i++){
        Matrix[i] = new double[2];
    }
    for(int i=0; i<n_datos;i++){
        for(int j=0;j<2;j++){
            Matrix[i][j]=0;
        }
    }
    for(int i=0; i<n_datos;i++){
        Matrix[i][0]=new_x[i];
    }
    for(int i=0; i<n_datos;i++){
        double pol=0;
        for(int j=0;j<n_datos;j++){
            double base_l=1;
            for(int l=0;l<n_datos;l++){
                if(l!=j){
                    base_l*=(((Matrix[i][0]-datos[l][0]))/(datos[j][0]-datos[l][0]));
                }
            }
            pol+=base_l*datos[j][1];
        }
        Matrix[i][1]=pol;
    }
    return Matrix;
}
double** TransformadaFourier(int n_datos,double** datos){
    double samplerate=(datos[0][0]-datos[n_datos-1][0])/(n_datos);
    double* freq=linspace(datos[0][0]*(samplerate/n_datos),datos[n_datos-1][0]*(samplerate/n_datos),n_datos);
    double* RTF = new double[n_datos];
    double* ITF = new double[n_datos];
    for(int i = 0; i < n_datos; i++)
    {
        complex<double> sum(0.0,0.0);
        for( int j = 0; j < n_datos; j++)
        {
            complex<double> arg(0.0, M_PI/n_datos*(-2)*j*i);
            sum += datos[j][1] *exp(arg);
        }
        RTF[i]= sum.real();
        ITF[i]= sum.imag();
    }
    double **Matrix = new double*[n_datos];
    for(int i=0;i<n_datos;i++){
        Matrix[i] = new double[3];
    }
    for(int i=0; i<n_datos;i++){
        Matrix[i][0]=freq[i];
        Matrix[i][1]=RTF[i];
        Matrix[i][2]=ITF[i];
    }
    return Matrix;
}
int main () {
    int NDatos;
    NDatos=Number_data();
    cout<< NDatos<<endl;
    double** Datos=DatatoMatrix(NDatos,2);
    double** interpolacion=PolinomiodeLagrange(NDatos,Datos);
    double** transform=TransformadaFourier(NDatos,interpolacion);
    ofstream fs("transformada.txt");
    for(int i=0; i<NDatos;i++){
        for(int j=0;j<3;j++){
            fs<<transform[i][j]<<" ";
        }
        fs<<" "<<endl;
    }
    fs.close();
    return 0;
}


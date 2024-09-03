#include<iostream>
#include<cmath>
#include<iomanip>
#include<fstream>
#include <stdio.h>
#include<stdlib.h>
using namespace std;

double Looped_Legendre_polynomial(double x, int N, double result[3]);
double Looped_Lobatto_polynomial(double x, int N, double output[2]);
double Newton_Raphson_method(double x, int N, double epsilon);

double Lagrange_Polynomials(int N, double *roots, double x, int i);
double first_derivative_Lagrange_Polynomials(int N, double *roots, double x, int i);

double Mass_Matrix_GQI(int N, int N_GQI, double *roots, double *weights, double *GQI_roots, double *GQI_weights, int k, int l );
double Laplacian_Matrix_GQI(int N, int N_GQI, double *roots, double *weights, double *GQI_roots, double *GQI_weights, int k, int l );
double Differentiation_Matrix_GQI(int N, int N_GQI, double *roots, double *weights, double *GQI_roots, double *GQI_weights, int k, int l );

double Gauss_Elimination(int p, double *H, double *rhs, double *q);

int main() {
    int e=4; // Number of elements
    int N=2; // number of nodal intervals in each element / degree of polynomial used for interpolating element

    fstream input_file;
    string t3;
    int t1;
    input_file.open("1D_Heat_Cond_input.dat",ios::in);
    // if (input_file == NULL) {
	// 	cout << "Can't Open input file !!!\n";
	// 	exit(1);
	// }
    input_file >> t1;
    e = t1;
    getline(input_file,t3);getline(input_file,t3);
    double element_size[e];
    cout << "e = " << e << endl;
    for(int i=0;i<e;i++) {
        input_file >> t1;
        cout << "t1 = " << t1 << endl;
        element_size[i] = t1;
        getline(input_file,t3);
        cout << "i = " << i << " -->  element_size = " << element_size[i] << endl;
    }
    input_file.close();

    /*
    FILE *fp1;
    fp1=fopen("1D_Heat_Cond.dat","r");
    if(fp1==NULL) {
        cout << "File is Empty or does not exist" << endl;
        exit(0);
    }
    fscanf(fp1,"%d",&e);
    double element_size[e];
    cout << "e = " << e << endl;
    for(int i=0;i<e;i++) {
        fscanf(fp1,"%lf",&element_size[i]);
        cout << "i = " << i << " -->  element_size = " << element_size[i] << endl;
    }
    fclose(fp1);*/

    /*double element_size[] = {2,2,1,4};
    for(int i=0;i<e;i++) {
        cout << "i = " << i << " -->  element_size = " << element_size[i] << endl;
    }*/


    // Generating the N+1 roots of the Lobatto Polynomial
    double roots[N+1];
    double weights[N+1];
    double q_start = 1;
    double q_end = 1;
    double temp; // temporary variable
    double epsilon = pow(10,-6);
    double result[3];
    for(int i=0;i<=N;i++) {
        double x = cos((double)(2*i)*M_PI/(2*N));
        // cout << "x = " << x << endl;
        roots[i]=Newton_Raphson_method(x, N+1, epsilon );
        // cout << roots[i] << endl;
        Looped_Legendre_polynomial(roots[i], N, result);
        temp=result[0];
        weights[i] = 2.0/(N*(N+1)*pow(temp,2));
    }
    cout << "The roots and weights of the Lobatto Polynomial for matrices are " << endl;
    for(int i=0;i<=N;i++) {
        cout << roots[i] << "\t" << weights[i] << endl;
    }
    cout << endl;

    int N_GQI = N+1; // For Exact Integration using GQI - Gaussian Quadrature Integration
    double GQI_roots[N_GQI+1];
    double GQI_weights[N_GQI+1];
    for(int i=0;i<=N_GQI;i++) {
        double x = cos((double)(2*i)*M_PI/(2*N_GQI));
        // cout << "x = " << x << endl;
        GQI_roots[i]=Newton_Raphson_method(x, N_GQI+1, epsilon );
        // cout << roots[i] << endl;
        Looped_Legendre_polynomial(GQI_roots[i], N_GQI, result);
        temp=result[0];
        GQI_weights[i] = 2.0/(N_GQI*(N_GQI+1)*pow(temp,2));
    }
    cout << "The roots and weights of the Lobatto Polynomial for GQI are " << endl;
    for(int i=0;i<=N_GQI;i++) {
        cout << GQI_roots[i] << "\t" << GQI_weights[i] << endl;
    }


    // Defining all Matrices
    int p = e*(N-1) + (e+1); // Total number of nodes
    double M[p][p]; // Mass Matrix
    double L[p][p]; // Laplacian Matrix
    double D[p][p]; // Differentiation Matrix
    // Initializing Mass and aplacian Matrix to zero
    for(int i=0;i<p;i++) {
        for(int j=0;j<p;j++) {
            M[i][j]=0;
            L[i][j]=0;
            D[i][j]=0;
            //cout << M[i][j] << "\t";
        }
        //cout << endl;
    }

    // Connectivity Matrix
    /*int connectivity_matrix[e][N+1];
    for(int i=0;i<e;i++) {
        for(int j=0;j<=N;j++) {
            connectivity_matrix[i][j]=i*e+j;
        }
    }
    // int I, J; // Global Coordinates
    // I=connectivity_matrix[g][k];
    // J=connectivity_matrix[g][l];
    */

    double Element_Matrix[N+1][N+1];
    int I, J; // Global Coordinates

    // Mass Matrix
    // cout << "Mass Matrix" << endl;
    for(int k=0;k<=N;k++) {
        for(int l=0;l<=N;l++) {
            Element_Matrix[k][l] = Mass_Matrix_GQI(N, N_GQI, roots, weights, GQI_roots, GQI_weights, k, l );
            // cout << Element_Matrix[k][l] << "\t"; // Testing
        }
        // cout << endl; // Testing
    }
    for(int g=0;g<e;g++) {
        for(int k=0;k<=N;k++) {
            for(int l=0;l<=N;l++) {
                I=(g)*(N) + k;
                J=(g)*(N) + l;
                temp = ((element_size[g])/2.0)*Element_Matrix[k][l];
                M[I][J] = M[I][J] + temp;
            }
        }
    }

    // Laplacian Matrix
    // cout << "Laplacian Matrix" << endl;
    for(int k=0;k<=N;k++) {
        for(int l=0;l<=N;l++) {
            Element_Matrix[k][l] = Laplacian_Matrix_GQI(N, N_GQI, roots, weights, GQI_roots, GQI_weights, k, l );
            // cout << Element_Matrix[k][l] << "\t"; // Testing
        }
        // cout << endl; // Testing
    }
    for(int g=0;g<e;g++) {
        for(int k=0;k<=N;k++) {
            for(int l=0;l<=N;l++) {
                I=(g)*(N) + k;
                J=(g)*(N) + l;
                temp = -(2.0/(element_size[g]))*Element_Matrix[k][l];
                L[I][J] = L[I][J] + temp;
            }
        }
    }

    // Differentiation Matrix
    // cout << "Differentiation Matrix" << endl;
    for(int k=0;k<=N;k++) {
        for(int l=0;l<=N;l++) {
            Element_Matrix[k][l] = Differentiation_Matrix_GQI(N, N_GQI, roots, weights, GQI_roots, GQI_weights, k, l );
            // cout << Element_Matrix[k][l] << "\t"; // Testing
        }
        // cout << endl; // Testing
    }
    for(int g=0;g<e;g++) {
        for(int k=0;k<=N;k++) {
            for(int l=0;l<=N;l++) {
                I=(g)*(N) + k;
                J=(g)*(N) + l;
                temp = (1.0)*Element_Matrix[k][l];
                D[I][J] = D[I][J] + temp;
            }
        }
    }

    // Displaying Matrix
    cout << "\n \n" <<  "**************** Mass Matrix *********************\n\n";
    for(int i=0;i<p;i++) {
        for(int j=0;j<p;j++) {
            if (M[i][j] == 0) {
                cout << "0.000000" << "  ";
            }
            else {
                cout << setprecision(6) << M[i][j] << "  ";
            }

        }
        cout << endl;
    }

    cout << "\n \n" << "**************** Laplacian Matrix ******************* \n\n";
    for(int i=0;i<p;i++) {
        for(int j=0;j<p;j++) {
            if (M[i][j] == 0) {
                cout << "0.000000" << "  ";
            }
            else {
                cout << setprecision(6) << L[i][j] << "  ";
            }
        }
        cout << endl;
    }

    cout << "\n \n" << "**************** Differentiation Matrix ******************* \n\n";
    for(int i=0;i<p;i++) {
        for(int j=0;j<p;j++) {
            if (D[i][j] == 0) {
                cout << "0.000000" << "  ";
            }
            else {
                cout << setprecision(6) << D[i][j] << "  ";
            }
        }
        cout << endl;
    }



    // RHS
    double f[p], x[p];
    double x_temp=0;
    x[0]=0;
    f[0]=0;
    cout << "x = " << x[0] << "\t f = " << f[0] << endl;
    int k=1;
    for(int g=0;g<e;g++) {
        for(int i=0;i<N;i++) {
            x_temp = x_temp + element_size[g]/N;
            x[k] = x_temp;
            f[k] = -pow(M_PI,2)*sin(M_PI*x[k]/9.0)/81.0;
            // cout << "x = " << x[k] << "\t f = " << f[k] << endl; // Testing
            k++;
        }
    }

    double rhs[p];
    // cout << "RHS" << endl;
    for(int i=0;i<p;i++) {
        rhs[i]=0;
        for(int j=0;j<p;j++) {
            rhs[i] = rhs[i] + M[i][j]*f[j];
        }
        // cout << rhs[i] << "\t";
    }
    // cout << endl;

    // Eliminating the Corner nodes from Matrix Operation
    for(int j=0;j<p;j++) {
        L[0][j] = 0;
        L[p-1][j] = 0;
    }
    for(int i=0;i<p;i++) {
        rhs[i] = rhs[i] - L[i][0]*q_start - L[i][p-1]*q_end;
        L[i][0] = 0;
        L[i][p-1] = 0;
    }
    L[0][0] = 1.0;
    L[p-1][p-1] = 1.0;
    rhs[0]=q_start;
    rhs[p-1]=q_end;

    // Testing
    /*cout << "Displaying Laplacian matrix after eliminating corner nodes" << endl;
    for(int i=0;i<p;i++) {
        for(int j=0;j<p;j++) {
            cout << L[i][j] << "  ";
        }
        cout << "\t | \t" << rhs[i] << endl;
        //cout << endl;
    }*/

    // Gauss Elimination
    // For passing 2d Array as 1D Array in C++
    double H[p*p];
    k=0;
    for(int i=0;i<p;i++) {
        for(int j=0;j<p;j++) {
            H[k] = L[i][j];
            k++;
        }
    }
    double q[p];
    Gauss_Elimination(p, H, rhs, q);

    cout << "x \t" << " Numerical_q \t" << " Actual_q" << endl;
    for(int i=0;i<p;i++) {
        cout << "x= " << x[i] << "\t" << "q_numerical = " << q[i] << " \t q_actual = " << sin(M_PI*x[i]/9.0) + 1.0 << endl;
    }

    //// ******** Printing Results ********
    fstream file;

    file.open("1D_Heat_Cond_Num_Simul.dat",ios::out);
    file << "x" << "\t" << "q_Num_simul"<< endl << endl;
    for(int i=0;i<p;i++) {
        file << x[i] << "\t" << q[i] << endl;;
    }
    file.close();

    file.open("1D_Heat_Cond_Actual_Result.dat",ios::out);
    file << "x" << "\t" << "q_actual"<< endl << endl;
    for(int i=0;i<p;i++) {
        file << x[i] << "\t" << sin(M_PI*x[i]/9.0) + 1.0 << endl;;
    }
    file.close();



    return 0;
}

double Gauss_Elimination(int p, double *H, double *rhs, double *q) {
    double A[p][p+1];
	for(int i=0;i<p;i++) {
		for(int j=0;j<p;j++) {
			// A[i][j] = L[i][j];
            A[i][j] = H[i*p + j];
		}
        A[i][p] = rhs[i];
    }

    // Testing
    /*for(int i=0;i<p;i++) {
        for(int j=0;j<p;j++) {
            cout << A[i][j] << "\t";
        }
        cout << "\t | \t " << A[i][p] << endl;
    }*/

	double temp, factor, sum, max;
	int swap_index;

	//Gauss Elimination
	for(int j=0;j<p-1;j++) {

		// Partial Pivoting
		max=fabs(A[j][j]);
        swap_index=j;
		for(int i=j+1;i<p;i++) {
			if (fabs(A[i][j])>max) {
				swap_index=i;
				max=fabs(A[i][j]);
			}
		}
		int l=swap_index;
		for(int i=j;i<=p;i++) {
			temp = A[j][i];
			A[j][i] = A[l][i];
			A[l][i] = temp;
		}

		// Forward Elimination
		if (max==0) {
			cout << "\n-------------The Matrix is Singular and doesn't have unique solution !!!!-------------\n" << endl;
			exit(0);
		}
		else {
			for(int i=j+1;i<p;i++) {
				factor=A[i][j]/A[j][j];
				for(int k=j;k<=p;k++) {
					A[i][k] = A[i][k] - factor*A[j][k];
				}
			}
		}

	}

	// Backward Substitution
	q[p-1]=A[p-1][p]/A[p-1][p-1];
	for(int i=p-2;i>=0;i--) {
		sum=0;
		for(int k=i+1;k<p;k++) {
			sum=sum + A[i][k]*q[k];
		}
		q[i]=(A[i][p] - sum)/A[i][i]; ;
	}

	return 1;
}

double Mass_Matrix_GQI(int N, int N_GQI, double *roots, double *weights, double *GQI_roots, double *GQI_weights,  int k, int l ) { // GQI = Gauss_Qudrature_integration
    double Integral=0;
    double f1, f2;
    for(int i=0;i<=N_GQI;i++) {
            f1=Lagrange_Polynomials(N, roots, GQI_roots[i], k);
            // cout << "f1 = " << f1 << endl;
            f2=Lagrange_Polynomials(N, roots, GQI_roots[i], l);
            // cout << "f2 = " << f2 << endl;
            Integral = Integral + GQI_weights[i]*f1*f2;
    }
    //cout << "Integral=" << Integral << endl;
    return Integral;
}

double Laplacian_Matrix_GQI(int N, int N_GQI, double *roots, double *weights, double *GQI_roots, double *GQI_weights,  int k, int l  ) { // GQI = Gauss_Qudrature_integration
    double Integral=0;
    double f1, f2;
    for(int i=0;i<=N_GQI;i++) {
        f1=first_derivative_Lagrange_Polynomials(N, roots, GQI_roots[i], k);
        f2=first_derivative_Lagrange_Polynomials(N, roots, GQI_roots[i], l);
        Integral = Integral + GQI_weights[i]*f1*f2;
    }
    return Integral;
}

double Differentiation_Matrix_GQI(int N, int N_GQI, double *roots, double *weights, double *GQI_roots, double *GQI_weights, int k, int l ) {
    double Integral=0;
    double f1, f2;
    for(int i=0;i<=N_GQI;i++) {
        f1=Lagrange_Polynomials(N, roots, GQI_roots[i], k);
        f2=first_derivative_Lagrange_Polynomials(N, roots, GQI_roots[i], l);
        Integral = Integral + GQI_weights[i]*f1*f2;
    }
    return Integral;
}

double Lagrange_Polynomials(int N, double *roots, double x, int s) {
    double L = 1;
    for(int j=0;j<=N;j++) {
        if(j!=s) {
            L = L*(x-roots[j])/(roots[s]-roots[j]);
        }
    }
    return L;
}

double first_derivative_Lagrange_Polynomials(int N, double *roots, double x, int i) {
    double L = 1;
    double S = 0;
    for(int k=0;k<=N;k++) {
        if (k!=i) {
            L=1;
            for(int j=0;j<=N;j++) {
                if(j!=i && j!=k) {
                    L = L*(x-roots[j])/(roots[i]-roots[j]);
                }
            }
            S = S + (1/(roots[i]-roots[k]))*L;
        }
    }
    return S;
}

double Looped_Legendre_polynomial(double x, int N, double result[3]) {
    double P = 0;
    double P_x = 0;
    double P_xx = 0;
    if (N==0) {
        P=1;
        P_x=0;
        P_xx=0;
    }
    else if (N==1) {
        P=x;
        P_x=1;
        P_xx=0;
    }
    if (N>=2) {
        // Polynomial value
        double a = 1;
        double b = x;
        // First Derivative value
        double c = 0;
        double d = 1;
        // Second Derivative value
        double e = 0;
        double f = 0;
        for(int i=2;i<=N;i++) {
            // Polynomial value
            P =  ((double)(2*i-1)/i)*x*b - ((double)(i-1)/i)*a;
            // First Derivative value
            P_x = ((double)(2*i-1)/i)*b +  ((double)(2*i-1)/i)*x*d - ((double)(i-1)/i)*c;
            // Second Derivative value
            P_xx = 2.0*((double)(2*i-1)/i)*d +  ((double)(2*i-1)/i)*x*f - ((double)(i-1)/i)*e;
            // ------value updates-----------
            // Polynomial value
            a = b;
            b = P;
            // First Derivative value
            c = d;
            d = P_x;
            // Second Derivative value
            e = f;
            f = P_xx;
        }
    }
    result[0]=P;
    result[1]=P_x;
    result[2]=P_xx;
    // cout << "Legendre Polynomial" << endl;
    // cout << P << "\t" << P_x << "\t" << P_xx << endl;

    return 0;
}

double Looped_Lobatto_polynomial(double x, int N, double output[2]) {
    double result[3];
    Looped_Legendre_polynomial(x, N-1, result);
    output[0]=(1.0-pow(x,2.0))*result[1];
    output[1]=(-2.0*x)*result[1] + (1-pow(x,2))*result[2];
    return 0;
}

double Newton_Raphson_method(double x, int N, double epsilon ) {
    double h=1;
    double output[2];
    while(fabs(h)>epsilon) {
        Looped_Lobatto_polynomial(x, N, output);
        h = output[0]/output[1];
        // cout << "h=" << h << endl;
        x = x - h;
    }
    // cout << x << endl;
    return x;
}

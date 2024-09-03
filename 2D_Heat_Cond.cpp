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

double Mass_Matrix_GQI(int Nx, int Ny, int N_GQI, double *roots_x, double *weights_x, double *roots_y, double *weights_y, double *GQI_roots, double *GQI_weights, int i, int j, int k, int l );
double Laplacian_Matrix_GQI(int Nx, int Ny, int N_GQI, double *roots_x, double *weights_x, double *roots_y, double *weights_y, double *GQI_roots, double *GQI_weights, int axis, int i, int j, int k, int l );
double Differentiation_Matrix_GQI( int Nx, int Ny, int N_GQI, double *roots_x, double *weights_x, double *roots_y, double *weights_y, double *GQI_roots, double *GQI_weights, int axis, int i, int j, int k, int l );

double Gauss_Elimination(int p, double *H, double *rhs, double *q);

int main() {
    int ex; // Number of elements in x direction
    int ey; // Number of elements in y direction
    int Nx=6; // number of nodal intervals in each element in x direction / order of polynomial used for interpolating element
    int Ny=6; // number of nodal intervals in each element in y direction / order of polynomial used for interpolating element
    
    // double element_size_x[] = {2,2,1,4};
    // double element_size_y[] = {2,2,1,4};

    double element_size_x[] = {2,2};
    ex = sizeof(element_size_x) / sizeof(element_size_x[0]); // Number of elements in x direction
    double element_size_y[] = {2,2};
    ey = sizeof(element_size_y) / sizeof(element_size_y[0]); // Number of elements in y direction
        
    // Generating the N+1 roots of the Lobatto Polynomial
    double roots_x[Nx+1];
    double weights_x[Nx+1];
    double roots_y[Ny+1];
    double weights_y[Ny+1];
    double q_top = 1.0;
    double q_bottom = 0.0;
    double q_right = 0.0;
    double q_left = 0.0;
    double temp, temp_1; // temporary variable
    double epsilon = pow(10,-6);
    double result[3];

    // For Polynomial Roots in x direction
    for(int i=0;i<=Nx;i++) {
        double x = cos((double)(2*i)*M_PI/(2*Nx));
        // cout << "x = " << x << endl; 
        roots_x[Nx-i]=Newton_Raphson_method(x, Nx+1, epsilon );
        // cout << roots_x[i] << endl;
        Looped_Legendre_polynomial(roots_x[Nx-i], Nx, result);
        temp=result[0];
        weights_x[i] = 2.0/((Nx)*(Nx+1)*pow(temp,2));
    }
    cout << "The roots and weights of the Lobatto Polynomial for FEM matrices in x direction are " << endl;
    for(int i=0;i<=Nx;i++) {
        cout << roots_x[i] << "\t" << weights_x[i] << endl;
    }
    cout << endl;

    // For Polynomial Roots in y direction
    for(int i=0;i<=Ny;i++) {
        double y = cos((double)(2*i)*M_PI/(2*Ny));
        // cout << "y = " << y << endl; 
        roots_y[Ny-i]=Newton_Raphson_method(y, Ny+1, epsilon );
        // cout << roots_y[i] << endl;
        Looped_Legendre_polynomial(roots_y[Ny-i], Ny, result);
        temp=result[0];
        weights_y[i] = 2.0/(Ny*(Ny+1)*pow(temp,2));
    }
    cout << "The roots and weights of the Lobatto Polynomial for FEM matrices in y direction are " << endl;
    for(int i=0;i<=Ny;i++) {
        cout << roots_y[i] << "\t" << weights_y[i] << endl;
    }
    cout << endl;

    // Roots/Points For Exact Integration using GQI - Gaussian Quadrature Integration
    int N_GQI = max(Nx,Ny)+1; 
    double GQI_roots[N_GQI+1];
    double GQI_weights[N_GQI+1];
    for(int i=0;i<=N_GQI;i++) {
        double x = cos((double)(2*i)*M_PI/(2*N_GQI));
        // cout << "x = " << x << endl; 
        GQI_roots[N_GQI-i]=Newton_Raphson_method(x, N_GQI+1, epsilon );
        // cout << roots[i] << endl;
        Looped_Legendre_polynomial(GQI_roots[N_GQI-i], N_GQI, result);
        temp=result[0];
        GQI_weights[i] = 2.0/(N_GQI*(N_GQI+1)*pow(temp,2));
    }
    cout << "The roots and weights of the Lobatto Polynomial for GQI are " << endl;
    for(int i=0;i<=N_GQI;i++) {
        cout << GQI_roots[i] << "\t" << GQI_weights[i] << endl;
    }


    // Defining all Matrices
    int p = ( ex*(Nx-1) + (ex+1) )*( ey*(Ny-1) + (ey+1) ); // Total number of nodes
    cout << "Total No. of nodes = " << p << endl;
    double M[p][p]; // Mass Matrix
    double L[p][p]; // Laplacian Matrix
    double D[p][p]; // Differentiaition Matrix
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

    // int I, J; // Global Coordinates
    
   // ex*Nx + 1 are total points in x direction
   // ey*Ny + 1 are total points in y direction
   // Connectivity Matrix
   double connect_matrix[ex*ey][(Nx+1)*(Ny+1)];
   for(int j=0;j<ey;j++) {
    for(int i=0;i<ex;i++) {
        int element = (j)*ex+i;
        for(int k=0;k<=Ny;k++) {
            for(int l=0;l<=Nx;l++) {
                int global_node = (j)*(ex*Nx+1)*Ny + k*(ex*Nx+1) + (i*Nx) + (l);
                connect_matrix[element][k*(Nx+1) + l] = global_node;
            }
        }
    }
   }
   // Testing
   cout << "The Connectivity Matrix for all elements is as follows" << endl;
   for(int j=0;j<ey;j++) {
    for(int i=0;i<ex;i++) {
        int element = (j)*ex+i;
        for(int k=0;k<=Ny;k++) {
                for(int l=0;l<=Nx;l++) {
                    cout << connect_matrix[element][k*(Nx+1) + l] << "\t";
                }
        }
        cout << endl;
    }
   }

    
    int EM_size = (Nx+1)*(Ny+1); // Element Matrix size
    double Element_Matrix_1[EM_size][EM_size];
    double Element_Matrix_2[EM_size][EM_size];
    int I, J; // Global Coordinates
    
    // Mass Matrix
    cout << "Mass Matrix" << endl;
    // Element Matrix
    for(int j=0;j<=Ny;j++) {
        for(int i=0;i<=Nx;i++) {
            for(int l=0;l<=Ny;l++) {
                for(int k=0;k<=Nx;k++) {
                    Element_Matrix_1[j*(Nx+1)+i][l*(Nx+1)+k] = Mass_Matrix_GQI(Nx, Ny, N_GQI, roots_x, weights_x, roots_y, weights_y, GQI_roots, GQI_weights, i, j, k, l );
                    cout << Element_Matrix_1[j*(Nx+1)+i][l*(Nx+1)+k] << "\t"; // Testing
                }
            }
            cout << endl; // Testing
        }
    }
    // Mass Global Matrix
    for(int h=0;h<ey;h++) {
        for(int g=0;g<ex;g++) {
            double Jacobian = element_size_x[g]*element_size_y[h]/4.0;
            int element = (h)*ex+g;
            for(int j=0;j<=Ny;j++) {
                for(int i=0;i<=Nx;i++) {
                    for(int l=0;l<=Ny;l++) {
                        for(int k=0;k<=Nx;k++) {
                            I=connect_matrix[element][j*(Nx+1)+i];
                            J=connect_matrix[element][l*(Nx+1)+k];
                            temp = Jacobian*Element_Matrix_1[j*(Nx+1)+i][l*(Nx+1)+k];
                            M[I][J] = M[I][J] + temp;
                        }
                    }
                }
            }
        }
    }

    
    // Laplacian Matrix
    cout << "Laplacian Matrix" << endl;
    // Element Matrix
    for(int j=0;j<=Ny;j++) {
        for(int i=0;i<=Nx;i++) {
            for(int l=0;l<=Ny;l++) {
                for(int k=0;k<=Nx;k++) {
                    Element_Matrix_1[j*(Nx+1)+i][l*(Nx+1)+k] = Laplacian_Matrix_GQI(Nx, Ny, N_GQI, roots_x, weights_x, roots_y, weights_y, GQI_roots, GQI_weights, 1, i, j, k, l );
                    // if (Element_Matrix_1[j*(Nx+1)+i][l*(Nx+1)+k]<pow(10,-3)) {     Element_Matrix_1[j*(Nx+1)+i][l*(Nx+1)+k] = 0.0;     }
                    Element_Matrix_2[j*(Nx+1)+i][l*(Nx+1)+k] = Laplacian_Matrix_GQI(Nx, Ny, N_GQI, roots_x, weights_x, roots_y, weights_y, GQI_roots, GQI_weights, 2, i, j, k, l );
                    // if (Element_Matrix_2[j*(Nx+1)+i][l*(Nx+1)+k]<pow(10,-3)) {     Element_Matrix_2[j*(Nx+1)+i][l*(Nx+1)+k] = 0.0;     }
                    // cout << Element_Matrix_1[j*(Nx+1)+i][l*(Nx+1)+k] << "\t"; // Testing
                    // cout << Element_Matrix_2[j*(Nx+1)+i][l*(Nx+1)+k] << "\t"; // Testing
                    cout << Element_Matrix_1[j*(Nx+1)+i][l*(Nx+1)+k] + Element_Matrix_2[j*(Nx+1)+i][l*(Nx+1)+k] << "\t"; // Testing
                }
            }
            cout << endl; // Testing
        }
    }
    // Laplacian Global Matrix
    for(int h=0;h<ey;h++) {
        for(int g=0;g<ex;g++) {
            double Jacobian_1 = -element_size_y[h]/element_size_x[g];
            double Jacobian_2 = -element_size_x[g]/element_size_y[h];
            int element = (h)*ex+g;
            for(int j=0;j<=Ny;j++) {
                for(int i=0;i<=Nx;i++) {
                    for(int l=0;l<=Ny;l++) {
                        for(int k=0;k<=Nx;k++) {
                            I=connect_matrix[element][j*(Nx+1)+i];
                            J=connect_matrix[element][l*(Nx+1)+k];
                            temp = Jacobian_1*Element_Matrix_1[j*(Nx+1)+i][l*(Nx+1)+k] + Jacobian_2*Element_Matrix_2[j*(Nx+1)+i][l*(Nx+1)+k];
                            L[I][J] = L[I][J] + temp;
                        }
                    }
                }
            }
        }
    }

    // Differentiation Matrix
    cout << "Differentiation Matrix" << endl;
    // Element Matrix
    for(int j=0;j<=Ny;j++) {
        for(int i=0;i<=Nx;i++) {
            for(int l=0;l<=Ny;l++) {
                for(int k=0;k<=Nx;k++) {
                    Element_Matrix_1[j*(Nx+1)+i][l*(Nx+1)+k] = Differentiation_Matrix_GQI(Nx, Ny, N_GQI, roots_x, weights_x, roots_y, weights_y, GQI_roots, GQI_weights, 1, i, j, k, l );
                    // if (Element_Matrix_1[j*(Nx+1)+i][l*(Nx+1)+k]<pow(10,-3)) {     Element_Matrix_1[j*(Nx+1)+i][l*(Nx+1)+k] = 0.0;     }
                    Element_Matrix_2[j*(Nx+1)+i][l*(Nx+1)+k] = Differentiation_Matrix_GQI(Nx, Ny, N_GQI, roots_x, weights_x, roots_y, weights_y, GQI_roots, GQI_weights, 2, i, j, k, l );
                    // if (Element_Matrix_2[j*(Nx+1)+i][l*(Nx+1)+k]<pow(10,-3)) {     Element_Matrix_2[j*(Nx+1)+i][l*(Nx+1)+k] = 0.0;     }
                    // cout << Element_Matrix_1[j*(Nx+1)+i][l*(Nx+1)+k] << "\t"; // Testing
                     cout << Element_Matrix_2[j*(Nx+1)+i][l*(Nx+1)+k] << "\t"; // Testing
                    // cout << Element_Matrix_1[j*(Nx+1)+i][l*(Nx+1)+k] + Element_Matrix_2[j*(Nx+1)+i][l*(Nx+1)+k] << "\t"; // Testing
                }
            }
            cout << endl; // Testing
        }
    }
    // Differentiation Global Matrix
    for(int h=0;h<ey;h++) {
        for(int g=0;g<ex;g++) {
            double Jacobian_1 = element_size_y[h]/2.0;
            double Jacobian_2 = element_size_x[g]/2.0;
            int element = (h)*ex+g;
            for(int j=0;j<=Ny;j++) {
                for(int i=0;i<=Nx;i++) {
                    for(int l=0;l<=Ny;l++) {
                        for(int k=0;k<=Nx;k++) {
                            I=connect_matrix[element][j*(Nx+1)+i];
                            J=connect_matrix[element][l*(Nx+1)+k];
                            temp = Jacobian_1*Element_Matrix_1[j*(Nx+1)+i][l*(Nx+1)+k] + Jacobian_2*Element_Matrix_2[j*(Nx+1)+i][l*(Nx+1)+k];
                            D[I][J] = D[I][J] + temp;
                        }
                    }
                }
            }
        }
    }
    
    /*
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
    */

    
    /////////////////// **************** Generating RHS Matrix and Position Coordinates ***********************
    double f[p], x[p], y[p];
    x[0]=0;
    y[0]=0;
    f[0]=0;
    int k=1;
    double N1, N2, zeta, x0,x1,y0,y1;
    for(int i=0;i<ex;i++) {
        x0=0;
        x1=element_size_x[i];
        temp=x[k-1];
        for(int j=1;j<=Nx;j++) {
            zeta=roots_x[j];
            N1 = (1-zeta)/2.0;
            N2 = (1+zeta)/2.0;
            x[k] = temp + N1*x0 + N2*x1;
            y[k] = 0;
            f[k]=0;
            k++;
        }
    }
    int x_nodes=ex*Nx+1; // Total nodes in x_direction
    int y_nodes=ey*Ny+1; // Total nodes in y_direction
    cout << endl << "x_nodes = " << x_nodes << "\t y_nodes = " << y_nodes << endl << endl;
    for(int i=0;i<ey;i++) {
        y0=0;
        y1=element_size_y[i];
        temp=y[k-1];
        for(int j=1;j<=Ny;j++) {
            zeta=roots_y[j];
            N1 = (1-zeta)/2.0;
            N2 = (1+zeta)/2.0;
            temp_1 = temp + N1*y0 + N2*y1;
            for(int l=0;l<x_nodes;l++) {
                x[k] = x[k-x_nodes];
                y[k] = temp_1;
                f[k]=0;
                k++;
            }
        }
    }
    // cout << "k = " << k << endl; // Testing
    // Testing
    //for(int i=0;i<p;i++) {
    //    cout << "x = " << x[i] << " y= " << y[i] << endl;
    //}
    
    
    
    //////// ************************* RHS of the System of Equations *******************
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

    ////////// *************** Eliminating the Boundary nodes from Matrix Operation ********************
    // Bottom Boundary
    for(int i=0;i<x_nodes;i++) { // For all entries in rows
        for(int j=0;j<p;j++) {
            L[i][j] = 0.0;
        }
        //L[i][i] = 1.0;
    }
    for(int i=0;i<p;i++) { // For all entries in columns
        for(int j=0;j<x_nodes;j++) {
            rhs[i] = rhs[i] - L[i][j]*q_bottom;
            // L[i][j]=0.0;
        }
    }
    // Left Boundary
    for(int i=0;i<p;i=i+x_nodes) { // For all entries in rows
        for(int j=0;j<p;j++) {
            L[i][j] = 0.0;
        }
        //L[i][i] = 1.0;
    }
    for(int i=0;i<p;i++) { // For all entries in columns
        for(int j=0;j<p;j=j+x_nodes) {
            rhs[i] = rhs[i] - L[i][j]*q_left;
            // L[i][j] = 0.0;
        }
    }
    // Right Boundary
    for(int i=x_nodes-1;i<p;i=i+x_nodes) { // For all entries in rows
        for(int j=0;j<p;j++) {
            L[i][j] = 0.0;
        }
        //L[i][i] = 1.0;
    }
    for(int i=0;i<p;i++) { // For all entries in columns
        for(int j=x_nodes-1;j<p;j=j+x_nodes) {
            rhs[i] = rhs[i] - L[i][j]*q_right;
            // L[i][j]=0.0;
        }
    }
    // Top Boundary
    for(int i=p-x_nodes;i<p;i++) { // For all entries in rows
        for(int j=0;j<p;j++) {
            L[i][j] = 0.0;
        }
        //L[i][i] = 1.0;
    }
    for(int i=0;i<p;i++) { // For all entries in columns
        for(int j=p-x_nodes;j<p;j++) {
            rhs[i] = rhs[i] - L[i][j]*q_top;
            // L[i][j]=0.0;
        }
    }
    ////////// To deal with Corner Nodes
    // Bottom Boundary
    for(int i=0;i<p;i++) { // For all entries in columns
        for(int j=0;j<x_nodes;j++) {
            L[i][j]=0.0;
        }
    }
    // Left Boundary
    for(int i=0;i<p;i++) { // For all entries in columns
        for(int j=0;j<p;j=j+x_nodes) {
            L[i][j] = 0.0;
        }
    }
    // Right Boundary
    for(int i=0;i<p;i++) { // For all entries in columns
        for(int j=x_nodes-1;j<p;j=j+x_nodes) {
            L[i][j]=0.0;
        }
    }
    // Top Boundary
    for(int i=0;i<p;i++) { // For all entries in columns
        for(int j=p-x_nodes;j<p;j++) {
            L[i][j]=0.0;
        }
    }



    ///////// *********** Setting RHS for all four boundaries as their respective values ***************
    // Bottom Boundary
    for(int i=0;i<x_nodes;i++) {
        L[i][i] = 1.0;
        rhs[i] = q_bottom;
    }
    // Left Boundary
    for(int i=0;i<p;i=i+x_nodes) {
        L[i][i] = 1.0;
        rhs[i] = q_left;
    }
    // Right Boundary
    for(int i=x_nodes-1;i<p;i=i+x_nodes) {
        L[i][i] = 1.0;
        rhs[i] = q_right;
    }
    // Top Boundary
    for(int i=p-x_nodes;i<p;i++) {
        L[i][i] = 1.0;
        rhs[i] = q_top;
    }

    

    // Testing
    cout << "Displaying Laplacian matrix after eliminating corner nodes" << endl;
    for(int i=0;i<p;i++) {
        for(int j=0;j<p;j++) {
            // if(L[i][j]<pow(10,-3)) {     L[i][j] = 0.0;      }
            cout << L[i][j] << "  ";
        }
        // if(rhs[i]<pow(10,-3)) {          rhs[i] = 0.0;       }
        cout << "\t | \t" << rhs[i] << endl;
        //cout << endl;
    }

    
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
        // cout << "x= " << x[i] << "\t" << "q_numerical = " << q[i] << " \t q_actual = " << sin(M_PI*x[i]/9.0) + 1.0 << endl;
        cout << "x= " << x[i] << "\t" << "q_numerical = " << q[i] <<  endl;
    }

    //// ******** Printing Results ********
    fstream file;

    file.open("2D_Heat_Cond_Num_Simul.dat",ios::out);
    file << "VARIABLES = \"X\",\"Y\",\"q_Num_simul\"\n";
    file << "ZONE T = \"BLOCK1\", I= " << x_nodes <<", J= " << y_nodes << ", F = POINT\n\n";
    for(int i=0;i<p;i++) {
        file << x[i] << "\t" << y[i] << "\t" << q[i] << endl;;
    }
    file.close();

    // file.open("2D_Heat_Cond_Actual_Result.dat",ios::out);
    // file << "VARIABLES =  = \"X\",\"Y\",\"q_actual\"\n";
    // file << "ZONE T = \"BLOCK1\", I= " << x_nodes <<", J= " << y_nodes << ", F = POINT\n\n";
    // for(int i=0;i<p;i++) {
    //     file << x[i] << y[i] <<"\t" << sin(M_PI*x[i]/9.0) + 1.0 << endl;;
    // }
    // file.close();
    
    

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

double Mass_Matrix_GQI(int Nx, int Ny, int N_GQI, double *roots_x, double *weights_x, double *roots_y, double *weights_y, double *GQI_roots, double *GQI_weights, int i, int j, int k, int l ) { // GQI = Gauss_Qudrature_integration
    double Integral=0;
    double f1, f2;
    //cout << i << j << k << l << endl;
    for(int r=0;r<=N_GQI;r++) {
        for(int s=0;s<=N_GQI;s++) {
            f1=Lagrange_Polynomials(Nx, roots_x, GQI_roots[r], i)*Lagrange_Polynomials(Ny, roots_y, GQI_roots[s], j);
            // cout << "f1 = " << f1 << endl;
            f2=Lagrange_Polynomials(Nx, roots_x, GQI_roots[r], k)*Lagrange_Polynomials(Ny, roots_y, GQI_roots[s], l);
            // cout << "f2 = " << f2 << endl;
            Integral = Integral + GQI_weights[r]*GQI_weights[s]*f1*f2;
            //cout << GQI_weights[r]*GQI_weights[s]*f1*f2 << endl;
        }
        //cout << endl;
    }
    //cout << "Integral=" << Integral << endl;
    return Integral;
}

double Laplacian_Matrix_GQI(int Nx, int Ny, int N_GQI, double *roots_x, double *weights_x, double *roots_y, double *weights_y, double *GQI_roots, double *GQI_weights, int axis, int i, int j, int k, int l ) { // GQI = Gauss_Qudrature_integration
    double Integral=0;
    double f1, f2;
    if (axis==1) {
        for(int r=0;r<=N_GQI;r++) {
            for(int s=0;s<=N_GQI;s++) {
                f1=first_derivative_Lagrange_Polynomials(Nx, roots_x, GQI_roots[r], i)*Lagrange_Polynomials(Ny, roots_y, GQI_roots[s], j);
                f2=first_derivative_Lagrange_Polynomials(Nx, roots_x, GQI_roots[r], k)*Lagrange_Polynomials(Ny, roots_y, GQI_roots[s], l);
                Integral = Integral + GQI_weights[r]*GQI_weights[s]*f1*f2;
            }
        }
    }
    if (axis==2) {
        for(int r=0;r<=N_GQI;r++) {
            for(int s=0;s<=N_GQI;s++) {
                f1=Lagrange_Polynomials(Nx, roots_x, GQI_roots[r], i)*first_derivative_Lagrange_Polynomials(Ny, roots_y, GQI_roots[s], j);
                f2=Lagrange_Polynomials(Nx, roots_x, GQI_roots[r], k)*first_derivative_Lagrange_Polynomials(Ny, roots_y, GQI_roots[s], l);
                Integral = Integral + GQI_weights[r]*GQI_weights[s]*f1*f2;
            }
        }
    }
    
    return Integral;
}

double Differentiation_Matrix_GQI(int Nx, int Ny, int N_GQI, double *roots_x, double *weights_x, double *roots_y, double *weights_y, double *GQI_roots, double *GQI_weights, int axis, int i, int j, int k, int l ) { // GQI = Gauss_Qudrature_integration
    double Integral=0;
    double f1, f2;
    if (axis==1) {
        for(int r=0;r<=N_GQI;r++) {
            for(int s=0;s<=N_GQI;s++) {
                f1=Lagrange_Polynomials(Nx, roots_x, GQI_roots[r], i)*Lagrange_Polynomials(Ny, roots_y, GQI_roots[s], j);
                f2=first_derivative_Lagrange_Polynomials(Nx, roots_x, GQI_roots[r], k)*Lagrange_Polynomials(Ny, roots_y, GQI_roots[s], l);
                Integral = Integral + GQI_weights[r]*GQI_weights[s]*f1*f2;
            }
        }
    }
    if (axis==2) {
        for(int r=0;r<=N_GQI;r++) {
            for(int s=0;s<=N_GQI;s++) {
                f1=Lagrange_Polynomials(Nx, roots_x, GQI_roots[r], i)*Lagrange_Polynomials(Ny, roots_y, GQI_roots[s], j);
                f2=Lagrange_Polynomials(Nx, roots_x, GQI_roots[r], k)*first_derivative_Lagrange_Polynomials(Ny, roots_y, GQI_roots[s], l);
                Integral = Integral + GQI_weights[r]*GQI_weights[s]*f1*f2;
            }
        }
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
            S = S + (1.0/(roots[i]-roots[k]))*L;
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
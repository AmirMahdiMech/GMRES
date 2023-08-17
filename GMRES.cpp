# include <cmath>
# include <cstdlib>
# include <fstream>
# include <iomanip>
# include <iostream>
# include <ctime>
using namespace std;

#define File "e20r0100.mtx"

void sizeread(int &size, int &NZ) {

	std::ifstream fin(File);
	
	// Declare variables:
	int M, N, L;
	
	// Ignore headers and comments:
	while (fin.peek() == '%') fin.ignore(2048, '\n');
	
	// Read defining parameters:
	fin >> M >> N >> L;
	size = M;
	// Create your matrix:
	double* matrix;		     // Creates a pointer to the array
	matrix = new double[M*N];	     // Creates the array of M*N size
	std::fill(matrix, matrix + M*N, 0.); // From <algorithm>, zeros all entries.
	
	// Read the data
	for (int l = 0; l < L; l++)
	{
		int m, n;
		double data;
		fin >> m >> n >> data;
		matrix[(m-1) + (n-1)*M] = data;
	}
	int counter = 0;
	for (int m = 0; m < M; m++)
	{
	    for(int n = 0; n < N; n++){
	        if (matrix[m + n*M]!=0)
	        counter++;
	        if (m == n && matrix[m + n*M]==0){
	        matrix[m + n*M] = 1E-6;
	        counter++;}
	}
	}
	NZ = counter;
	fin.close();	
}

void Matrixread(double a[], int ra[], int ca[]) {

	std::ifstream fin(File);
	
	int M, N, L;
	
	while (fin.peek() == '%') fin.ignore(2048, '\n');
	
	fin >> M >> N >> L;
	
	double* matrix;		     
	matrix = new double[M*N];	     
	std::fill(matrix, matrix + M*N, 0.); 
	
	for (int l = 0; l < L; l++)
	{
		int m, n;
		double data;
		fin >> m >> n >> data;
		matrix[(m-1) + (n-1)*M] = data;
	}
	int counter = 0;
	for (int m = 0; m < M; m++)
	{
	    for(int n = 0; n < N; n++){
	        if (matrix[m + n*M]!=0)
	        counter++;
	        if (m == n && matrix[m + n*M]==0){
	        matrix[m + n*M] = 1E-6;
	        counter++;}
	}
	}
	int c = 0;
    double threshold = 1E-7;
    //int max_nonzeros = M * N;
    //a =  new double[counter];
    //ra = new int[M+1];
	//ca = new int[counter];
    ra[0] = 0;
    for (int m = 0; m < M; m++) {
        for (int n = 0; n < N; n++) {
            double val = matrix[m * N + n];
            if (val >= threshold && c < counter) {
                a[c] = val;
                ca[c] = n;
                c++;
            }
        }
        ra[m+1] = c;
	}
	fin.close();
	return;	
}


double dotproduct ( int n, double a1[], double a2[] )
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a1[i] * a2[i];
  }
  return value;
}

void sort ( int n, int nz_num, int ia[], int ja[], double a[] )
{
  
  int i;
  int is;
  int dummy1;
  int j;
  int j1;
  int j2;
  int k;
  double dummy2;

  for ( i = 0; i < n; i++ )
  {
    j1 = ia[i];
    j2 = ia[i+1];
    is = j2 - j1;

    for ( k = 1; k < is; k++ ) 
    {
      for ( j = j1; j < j2 - k; j++ ) 
      {
        if ( ja[j+1] < ja[j] ) 
        {
          dummy1 = ja[j+1];
          ja[j+1] =  ja[j];
          ja[j] =  dummy1;

          dummy2 = a[j+1];
          a[j+1] =  a[j];
          a[j] = dummy2;
        }
      }
    }
  }
  return;
}  

void GivensRotation ( double c, double s, int k, double g[] )
{
  double g1;
  double g2;

  g1 = c * g[k] - s * g[k+1];
  g2 = s * g[k] + c * g[k+1];

  g[k]   = g1;
  g[k+1] = g2;

  return;
}

void SparseAx ( int n, int nz_num, int ia[], int ja[], double a[], double x[],
  double w[] )
{
  int i;
  int k;
  int k1;
  int k2;

  for ( i = 0; i < n; i++ )
  {
    w[i] = 0.0;
    k1 = ia[i];
    k2 = ia[i+1];
    for ( k = k1; k < k2; k++ )
    {
      w[i] = w[i] + a[k] * x[ja[k]];
    }
  }
  return;
}

void diagonalIndex ( int n, int nz_num, int ia[], int ja[], int Index[] )
{
  int i;
  int j;
  int j1;
  int j2;

  for ( i = 0; i < n; i++ )
  {
    Index[i] = -1;
    j1 = ia[i];
    j2 = ia[i+1];

    for ( j = j1; j < j2; j++ )
    {
      if ( ja[j] == i ) 
      {
        Index[i] = j;
      }
    }

  }
  return;
}

void ILU_factorization ( int n, int nz_num, int ia[], int ja[], double a[], int ua[],
  double l[] )
{
  int *iw;
  int i;
  int j;
  int jj;
  int jrow;
  int jw;
  int k;
  double tl;

  iw = new int[n];

  for ( k = 0; k < nz_num; k++ ) 
  {
    l[k] = a[k];
  }
  for ( i = 0; i < n; i++ ) 
  {

    for ( j = 0; j < n; j++ )
    {
      iw[j] = -1;
    }

    for ( k = ia[i]; k <= ia[i+1] - 1; k++ ) 
    {
      iw[ja[k]] = k;
    }

    j = ia[i];
    do 
    {
      jrow = ja[j];
      if ( i <= jrow )
      {
        break;
      }
      tl = l[j] * l[ua[jrow]];
      l[j] = tl;
      for ( jj = ua[jrow] + 1; jj <= ia[jrow+1] - 1; jj++ ) 
      {
        jw = iw[ja[jj]];
        if ( jw != -1 ) 
        {
          l[jw] = l[jw] - tl * l[jj];
        }
      }
      j = j + 1;
    } while ( j <= ia[i+1] - 1 );

    ua[i] = j;
    if ( jrow != i || l[j] == 0.0) 
    {
      cerr << "\n";
      cerr << "  Zero pivot on row " << i << "\n";
      exit ( 1 );
    }


    l[j] = 1.0 / l[j];
  }

  for ( k = 0; k < n; k++ ) 
  {
    l[ua[k]] = 1.0 / l[ua[k]];
  }
  delete [] iw;

  return;
}

void Preconditioner ( int n, int nz_num, int ia[], int ja[], double l[], int ua[], 
  double r[], double z[] )
{
  int i;
  int j;
  double *w;

  w = new double[n];

  for ( i = 0; i < n; i++ )
  {
    w[i] = r[i];
  }

  for ( i = 1; i < n; i++ )
  {
    for ( j = ia[i]; j < ua[i]; j++ )
    {
      w[i] = w[i] - l[j] * w[ja[j]];
    }
  }

  for ( i = n - 1; 0 <= i; i-- ) 
  {
    for ( j = ua[i] + 1; j < ia[i+1]; j++ ) 
    {
      w[i] = w[i] - l[j] * w[ja[j]];
    }
    w[i] = w[i] / l[ua[i]];
  }

  for ( i = 0; i < n; i++ )
  {
    z[i] = w[i];
  }

  delete [] w;

  return;
}


void PreconditionedGMRES ( int n, int nz_num, int ia[], int ja[], double a[], 
  double x[], double rhs[])
{
  double tol_rel = 1.0e-8;
  double tol_abs;
  double av;
  int itr_max;
  double *c;
  double delta = 1.0e-03;
  double *g;
  double *h;
  double dummy;
  int i;
  int itr;
  int itr_used=0;
  int j;
  int k;
  int k_copy;
  double *l;
  int mr;
  double mu;
  double *r;
  double rho;
  double rho_tol;
  bool ok = false;
  double *s;
  int *ua;
  double *v;
  double *y;
  cout<<" Enter the maximum outer iterations: ";
  cin>>itr_max;
  while (!ok){
  cout<<" Enter the maximum inner iteration(size of the krylov subspace): ";
  cin>>mr;
  if (mr > n ){
  cout<<" It must be less than the size of the matrix("<<n<<")"<<" Try again "<<endl;
  }
  else {
  	ok = true;
  }
	}
  cout<<" Enter your desired tolerance for the solution: ";
  cin>>tol_abs;
  c = new double[mr+1];
  g = new double[mr+1];
  h = new double[(mr+1)*mr];
//l = new double[ia[n]+1];
  l = new double[nz_num];
  r = new double[n];
  s = new double[mr+1];
  ua = new int[n];
  v = new double[(mr+1)*n];
  y = new double[mr+1];
  
  

  sort ( n, nz_num, ia, ja, a );

  diagonalIndex ( n, nz_num, ia, ja, ua );

  //ILU_factorization ( n, nz_num, ia, ja, a, ua, l );

  for ( itr = 0; itr < itr_max; itr++ ) 
  {
    SparseAx ( n, nz_num, ia, ja, a, x, r );

    for ( i = 0; i < n; i++ ) 
    {
      r[i] = rhs[i] - r[i];
    }

    //Preconditioner ( n, nz_num, ia, ja, l, ua, r, r );

    rho = sqrt ( dotproduct ( n, r, r ) );

    if ( itr == 0 )
    {
      rho_tol = rho * tol_rel;
    }

    for ( i = 0; i < n; i++ ) 
    {
      v[i] = r[i] / rho;
    }

    g[0] = rho;
    for ( i = 1; i < mr + 1; i++ ) 
    {
      g[i] = 0.0;
    }

    for ( i = 0; i < mr + 1; i++ ) 
    {
      for ( j = 0; j < mr; j++ ) 
      {
        h[i*(mr)+j] = 0.0;
      }
    }

    for ( k = 0; k < mr; k++ )
    {
      k_copy = k;

      SparseAx ( n, nz_num, ia, ja, a, v+k*n, v+(k+1)*n ); 

      //Preconditioner ( n, nz_num, ia, ja, l, ua, v+(k+1)*n, v+(k+1)*n );

      av = sqrt ( dotproduct ( n, v+(k+1)*n, v+(k+1)*n ) );

      for ( j = 0; j <= k; j++ ) 
      {
        h[j*mr+k] = dotproduct ( n, v+(k+1)*n, v+j*n );
        for ( i = 0; i < n; i++ ) 
        {
          v[(k+1)*n+i] = v[(k+1)*n+i] - h[j*mr+k] * v[j*n+i];
        }
      }
      h[(k+1)*mr+k] = sqrt ( dotproduct ( n, v+(k+1)*n, v+(k+1)*n ) );

      if ( ( av + delta * h[(k+1)*mr+k]) == av ) 
      {
        for ( j = 0; j < k + 1; j++ )
        {
          dummy = dotproduct ( n, v+(k+1)*n, v+j*n );
          h[j*mr+k] = h[j*mr+k] + dummy;
          for ( i = 0; i < n; i++ ) 
          {
            v[(k+1)*n+i] = v[(k+1)*n+i] - dummy * v[j*n+i];
          }
        }
        h[(k+1)*mr+k] = sqrt ( dotproduct ( n, v+(k+1)*n, v+(k+1)*n ) );
      }

      if ( h[(k+1)*mr+k] != 0.0 )
      {
        for ( i = 0; i < n; i++ )
        {
          v[(k+1)*n+i] = v[(k+1)*n+i] / h[(k+1)*mr+k];
        } 
      }

      if ( 0 < k )  
      {
        for ( i = 0; i < k + 2; i++ ) 
        {
          y[i] = h[i*mr+k];
        }
        for ( j = 0; j < k; j++ ) 
        {
          GivensRotation ( c[j], s[j], j, y );
        }
        for ( i = 0; i < k + 2; i++ )
        {
          h[i*mr+k] = y[i];
        }
      }
      mu = sqrt ( h[k*mr+k] * h[k*mr+k] + h[(k+1)*mr+k] * h[(k+1)*mr+k] );
      c[k] = h[k*mr+k] / mu;
      s[k] = -h[(k+1)*mr+k] / mu;
      h[k*mr+k] = c[k] * h[k*mr+k] - s[k] * h[(k+1)*mr+k];
      h[(k+1)*mr+k] = 0.0;
      GivensRotation ( c[k], s[k], k, g );

      rho = fabs ( g[k+1] );

      itr_used = itr_used + 1;

      if ( rho <= rho_tol && rho <= tol_abs )
      {
        break;
      }
    }

    k = k_copy;

    y[k] = g[k] / h[k*mr+k];
    for ( i = k - 1; 0 <= i; i-- )
    {
      y[i] = g[i];
      for ( j = i + 1; j < k + 1; j++ ) 
      {
        y[i] = y[i] - h[i*mr+j] * y[j];
      }
      y[i] = y[i] / h[i*mr+i];
    }
    for ( i = 0; i < n; i++ )
    {
      for ( j = 0; j < k + 1; j++ )
      {
        x[i] = x[i] + v[j*n+i] * y[j];
      }
    }

    if ( rho <= rho_tol && rho <= tol_abs )
    {
      break;
    }
    /*double error = 0.0;
    for (int i = 0; i < n; i++) {
    	error += (x[i]-1)*(x[i]-1);
    }*/
    //error = sqrt(error);
    cout<<"Iteration: "<<itr+1<<"\t"<<"Residual: "<<rho<<endl;
  }

    cout << "\n";;
    cout << "PreconditionedGMRES:\n";
    cout << "  Iterations = " << itr+1 << "\n";
    cout << "  Final residual = " << rho << "\n";
  
  delete [] c;
  delete [] g;
  delete [] h;
  delete [] l;
  delete [] r;
  delete [] s;
  delete [] ua;
  delete [] v;
  delete [] y;

  return;
}

void MTXtest() {
// This is the actual test on the matrix downloaded from Matrixmarket
  int N , NZ_NUM;
  sizeread(N,NZ_NUM);
  double a[NZ_NUM];
  int ra[N+1];
  int ca[NZ_NUM];
  Matrixread(a, ra, ca);
  int itr_max;
  int i;
  int mr;
  int n = N;
  int nz_num = NZ_NUM;
  double rhs[N];
  double x[N];
  double x_estimate[N];
  srand(time(NULL));
  for ( i = 0; i < n; i++ )
    {
      x[i] = 1.0;
      x_estimate [i]= rand()*1.0/RAND_MAX;
      //cout<<x_estimate[i]<<endl;
    }
  	SparseAx ( N, NZ_NUM, ra, ca, a, x, rhs);
  	int test;
  	double tol_abs;
  	double tol_rel;
  	double x_error;

  
  	cout << "\n";
  	cout << "  Testing Preconditioned_GMRES on a sample matrix from matrixmarket.\n";
  	cout << "\n";
    cout << "\n";
    cout << "  Matrix size is N = " << n << "\n";

    PreconditionedGMRES ( n, nz_num, ra, ca, a, x_estimate, rhs);
	cout<<"  Final Solution: "<<endl;
    for ( i = 0; i < n; i++ )
    {
      cout<<"\t x["<<i+1<<"] = "<<x[i]<<endl;
    }
  return;

} 
void matrixToCSR(int &n, int &counter, double** matrix, int *ra, int *ca, double *a){
	
	counter = 0;
	cout<<"Your matrix is "<<n<<" by "<<n<<endl;
	cout<<"Now you have to input the elements: "<<endl;
	
    matrix = new double*[n];
    bool ok = false;

    while (!ok){
    for (int i = 0; i < n; i++) {
        matrix[i] = new double[n];
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout<<"a_"<<i+1<<j+1<<": ";
            cin >> matrix[i][j];
            cout<<endl;
            if(matrix[i][j]!=0){
            	counter++;
			}
			else if(i == j){
				matrix[i][j] = 1E-6;
			}
        }
    }
    cout<<"Your matrix is :"<<endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout<<matrix[i][j]<<"\t";
            
        }
            cout<<endl;
    }
    cout<<"If this is your matrix press 1, if not press 0 :"<<"\t";
    cin>>ok;
	}
    ra[0] = 0;
    int k = 0;
    for (int i = 0; i < n; i++) {

        ra[i] = k;
        for (int j = 0; j < n; j++) {
            if (matrix[i][j] != 0.0) {
                ca[k] = j;
                a[k] = matrix[i][j];
                k++;
            }
        }
    }
    ra[n] = k;  
    delete matrix;
}
void Yourtest(){
	int counter;
    int n;
	cout<<" Please enter the size of your matrix: ";
	cin>>n;
    double **matrix;
    int Max = n*n;
    int ra[n+1];
    double a[Max];
    int ca[Max];
    double b[n];
    double x[n];
    matrixToCSR(n,counter,matrix,ra,ca,a);
    int nz_num = counter;
	cout<<" Input the RHS of the equation: "<<endl;
    for ( int i = 0; i < n; i++ )
    {
      cout<<"\t b["<<i+1<<"] = ";
      cin>>b[i];
      x[i] = 0.0;
    }
    PreconditionedGMRES ( n, nz_num, ra, ca, a, x, b);
    cout<<"  Final Solution: "<<endl;
    for ( int i = 0; i < n; i++ )
    {
      cout<<"\t x["<<i+1<<"] = "<<x[i]<<endl;
    }
}



int main(int argc, char** argv) {
	//Yourtest();
	MTXtest();
	return 0;
}




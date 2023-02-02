		/*This program is created to handle some basic
		  linear algebra operations based on 'Introduction 
		  To Linear Algebra' by professor Gilbert Strang. 
		  Hope this will be useful.*/
		  
		/*author: Benchang*/
#include<iostream>
#include<vector>
#include<math.h>
#include<iomanip>
#include<set>
using namespace std;
typedef vector< vector<double> > two_dvector;
/*functions*/

bool is_in(const vector<int>& v,int ele);/*get free columns*/

int sgn(long double);
void givens_rotation(double,double,long double&,long double&);

two_dvector resize(int r,int c);

void print_manu(ostream& outs);

class matrix{
	private:
		two_dvector co;/*entry of the matrix*/
		int row;
		int col;
		two_dvector low;/*L*/
		two_dvector up;/*U*/
		two_dvector d;/*D*/
		two_dvector permu;/*P*/
		two_dvector inver;
		int per=0;/*times of permutations when doing LU*/
		two_dvector rref;/*row reduced echelon form*/
		int r=0;/*rank*/
		double determine;
		
		two_dvector GQ;/*Givens Q*/
		two_dvector GR;/*Givens R*/
		
		two_dvector GSQ;/*Gram-Schmidt Q*/
		two_dvector GSR;/*Gram-Schmidt R*/
		
		vector<int> pivot_col;/*pivot_col from rref*/
		vector<int> free_col;/*free_col from rref*/
		
		/*four fundamental subspaces*/
		two_dvector null_space;
		two_dvector col_space;
		two_dvector row_space;
		two_dvector left_null_space;
		
		set<double, greater<double> > evalue;/*eigenvalue(no repitition)*/
		vector<double> evaluev;/*eigenvalue(with possible repitition)*/
		two_dvector evector;/*eigenvector*/
		
		two_dvector Q;/*orthonormal matrix from Q ¢GN QT*/
		two_dvector lambda;/*diagonal matrix from Q ¢GN QT*/
		two_dvector S;/*S ¢GN S-1*/
		
		two_dvector U;/*U ¢GUVT in column space*/
		two_dvector sigma;/*singular value matrix*/
		two_dvector V;/*U ¢GUVT in row space*/
		
		/*helper function*/
		void sort_for_svd();
	public:
		/*constructor*/
		matrix();
		matrix(int r,int c);
		matrix(two_dvector& v);
		/*assignment operator*/
		void operator=(const matrix& a);
		
		/*functions for input and output*/
		friend istream& operator>>(istream& ins,matrix& a);
		friend ostream& operator<<(ostream& outs,const matrix& a);
		friend ostream& operator<<(ostream& outs,const two_dvector& a);
		friend ostream& operator<<(ostream& outs,const set<double,greater<double>> &s);
		void getlu(ostream& outs);
		void getldu(ostream& outs);
		void getGQR(ostream& outs);
		void getGSQR(ostream& outs);
		void getQ_lambda_QT(ostream& outs);
		void getS_lambda_S1(ostream& outs);
		void getsvd(ostream& outs);
		
		/*operator for(+, -, *)*/
		/*matrix,matrix*/
		friend matrix operator+(const matrix& a,const matrix& b);
		friend matrix operator-(const matrix& a,const matrix& b);
		friend matrix operator*(const matrix& a,const matrix& b);
		
		/*two_dvector,two_dvector*/
		friend two_dvector operator+(const two_dvector&,const two_dvector&);
		friend two_dvector operator-(const two_dvector&,const two_dvector&);
		friend two_dvector operator*(const two_dvector&,const two_dvector&);
		
		
		/*vector,vector*/
		friend vector<double> operator+(const vector<double>&,const vector<double>&);
		friend vector<double> operator-(const vector<double>&,const vector<double>&);
		friend double operator*(const vector<double>&,const vector<double>&);
		
		/*const*/
		friend matrix operator*(const double&,const matrix&);
		friend vector<double> operator*(const double&,const vector<double>&);
		friend vector<double> operator/(const vector<double>&,const double&);
		
		/*two_dvector,vector*/
		friend vector<double> operator*(const two_dvector&,const vector<double>&);
		
		/*functions for two_dvector*/
		
		friend void set_to_identity(two_dvector& v);
		friend two_dvector trans(const two_dvector &v);
		friend void normalize_col(two_dvector& v);
		friend void orthonormalize_col(two_dvector& v);
		
		friend two_dvector inversev(const two_dvector& v);
		bool is_square();
		bool is_symmetric();
		
		/*functions for matrix*/
		int get_row();
		int get_col();
		
		void lu();
		void ldu();
		matrix transpose() const;
		matrix inverse();
		matrix get_rref();
		int rank();
		matrix project();
		double det();
		/*subspaces and orthogonolization*/
		two_dvector nullspace();
		two_dvector rowspace();
		two_dvector colspace();
		two_dvector leftnullspace();
		void GQR();
		void GSQR();
		
		/*eigenvalue and diagonalization*/
		set<double, greater<double> > eigenvalue();
		two_dvector eigenvector();
		void Q_lambda_QT();
		void S_lambda_S1();
		void svd();
		/*approximation*/
		matrix rank_n_appr(int n);
		/*Pseudo inverse*/
		matrix pseudo_inverse();
		/*last fundamental problem*/
		friend matrix Axb(matrix& a,matrix& b);
		/* condition number */
		double con_num();
		/*destructor*/
		~matrix();
};
////////////////////*constructors*//////////////////
matrix::matrix(){
	row=0;
	col=0;
}

matrix::matrix(int r,int c){
	row=r;
	col=c;
	co.resize(row);
	for(int i=0;i<co.size();i++)
		co[i].resize(col,0);
}

matrix::matrix(two_dvector& v){
	co=v;
	row=v.size();
	col=v[0].size();
}
//////////////*functions for input and output*///////////////

istream& operator>>(istream& ins,matrix& a){
	cout<<"Input rows and columns of the matrix"<<endl;
	ins>>a.row>>a.col;
	a.co.resize(a.row);
	for(int i=0;i<a.co.size();i++)
		a.co[i].resize(a.col,0);
	cout<<"Please input a "<<a.row<<" by "<<a.col<<" matrix"<<endl;
	for(auto& r: a.co)
		for(auto& c:r)
			ins>>c;
	return ins;
}

ostream& operator<<(ostream& outs,const two_dvector& v){
	for(auto r : v){
		for(auto c : r)
			if(fabs(c)<0.0001)
				outs<<right<<setw(10)<<"0"<<" ";
			else
				outs<<right<<setw(10)<<c<<" ";
		cout<<endl;
	}
	return outs;	
}

ostream& operator<<(ostream& outs,const set<double,greater<double>> &s){
	for(auto x:s){
		if(fabs(x)<0.0001)
			outs<<"0"<<" ";
		else
			outs<<x<<" ";
	}
	cout<<endl;
	return outs;
}

ostream& operator<<(ostream& outs,const matrix& a){
	for(auto r : a.co){
		for(auto c : r)
			if(fabs(c)<0.0001)
				outs<<right<<setw(10)<<"0"<<" ";
			else
				outs<<right<<setw(10)<<c<<" ";
		cout<<endl;
	}
	return outs;	
}

///////////////////*assignment operator*////////////////

void matrix::operator=(const matrix& a){
	row=a.row;
	col=a.col;
	per=a.per;
	r=a.r;
	determine=a.determine;
	co=a.co;
	low=a.low;
	up=a.up;
	d=a.d;
	rref=a.rref;
	permu=a.permu;
	inver=a.inver;
	null_space=a.null_space;
	col_space=a.col_space;
	row_space=a.row_space;
	GQ=a.GQ;
	GR=a.GR;
	GSQ=a.GSQ;
	GSR=a.GSR;
	evalue=a.evalue;
	evector=a.evector;
}
///////////////////*operator for(+, -, *)*//////////////

/*two_dvector,two_dvector*/

two_dvector operator+(const two_dvector& a,const two_dvector& b){
	two_dvector ans=resize(a.size(),a[0].size());
	for(int i=0;i<a.size();i++)
		for(int j=0;j<a[0].size();j++)
			ans[i][j]=(a[i][j]+b[i][j]);
	return ans;
}

two_dvector operator-(const two_dvector& a,const two_dvector& b){
	two_dvector ans=resize(a.size(),a[0].size());
	for(int i=0;i<a.size();i++)
		for(int j=0;j<a[0].size();j++)
			ans[i][j]=a[i][j]-b[i][j];
	return ans;
}

two_dvector operator*(const two_dvector& a,const two_dvector& b){
	two_dvector ans=resize(a.size(),b[0].size());
		for(int i=0;i<a.size();i++)
			for(int j=0;j<b[0].size();j++)
				for(int k=0;k<b.size();k++)
					ans[i][j]+=a[i][k]*b[k][j];
	return ans;
}


/*matrix,matrix*/

matrix operator+(const matrix& a,const matrix& b){
	if(a.row==b.row&&a.col==b.col){
		matrix ans(a.row,a.col);
		ans.co=a.co+b.co;
		return ans;
	}
	else{
		cout<<"Invalid input(rows and columns of two matrices should be the same)"<<endl;
		matrix ans(0,0);
		return ans;
	}
}
matrix operator-(const matrix& a,const matrix& b){
	if(a.row==b.row&&a.col==b.col){
		matrix ans(a.row,a.col);
		ans.co=a.co-b.co;
		return ans;
	}
	else{
		cout<<"Invalid input(rows and columns of two matrices should be the same)"<<endl;
		matrix ans(0,0);
		return ans;
	}
}

matrix operator*(const matrix& a,const matrix& b){
	if(a.col==b.row){
		matrix ans(a.row,b.col);
		ans.co=a.co*b.co;
		return ans;
	}
	else{
		cout<<"Invalid input(A's col should equal to B's row)"<<endl;
		matrix ans(0,0);
		return ans;
	}
}

/*vector,vector*/

vector<double> operator+(const vector<double>& a,const vector<double>& b){
	vector<double> ans(a.size(),0);
	if(a.size()!=b.size()){
		cout<<"Warning!!!vector addition of vector invalid"<<endl;
		return ans;
	}
	for(int i=0;i<a.size();i++)
		ans[i]+=a[i]+b[i];
	return ans;
}

vector<double> operator-(const vector<double>& a,const vector<double>& b){
	vector<double> ans(a.size(),0);
	if(a.size()!=b.size()){
		cout<<"Warning!!!vector substraction of vector invalid"<<endl;
		return ans;
	}
	for(int i=0;i<a.size();i++)
		ans[i]+=a[i]-b[i];
	return ans;
}
double operator*(const vector<double>& a,const vector<double>& b){
	double ans=0;
	if(a.size()!=b.size()){
		cout<<"Warning!!!inner product of vector invalid"<<endl;
		return 0;
	}
	for(int i=0;i<a.size();i++)
		ans+=a[i]*b[i];
	return ans;
}


/*const*/
matrix operator*(const double& n,const matrix& b){
	matrix ans=b;
	for(auto& x:ans.co)
		for(auto& y:x)
			y*=n;
	return ans;
}

two_dvector operator*(const double& n,const two_dvector &v){
	two_dvector ans=v;
	for(auto& x:ans)
		for(auto& y:x)
			y*=n;
	return ans;
}

vector<double> operator*(const double& n,const vector<double>& v){
	vector<double> ans=v;
	for(auto& x:ans)
		x*=n;
	return ans;
}
vector<double> operator/(const vector<double>& v,const double &n){
	vector<double> ans=v;
	for(auto& x:ans)
		x/=n;
	return ans;
}
/*two_dvector,vector*/

vector<double> operator*(const two_dvector& a,const vector<double>& v){
	vector<double> ans(a.size(),0);
	for(int i=0;i<a.size();i++)
		for(int j=0;j<v.size();j++)
			ans[i]+=a[i][j]*v[j];
	return ans;
}

//////////////*functions for two_dvector*////////////////
		
two_dvector resize(int r,int c){
	two_dvector ans;
	ans.resize(r);
	for(int i=0;i<ans.size();i++)
		ans[i].resize(c,0);
	return ans;
}

void set_to_identity(two_dvector& v){
	for(int i=0;i<v.size();i++)
		for(int j=0;j<v[0].size();j++)
			v[i][j]=(i==j);
}

two_dvector trans(const two_dvector &v){
	two_dvector ans;
	ans=resize(v[0].size(),v.size());
	for(int i=0;i<v.size();i++)
		for(int j=0;j<v[0].size();j++)
			ans[j][i]=v[i][j];
	return ans;	
}
void normalize_col(two_dvector& v){
	v=trans(v);
	double length;
	for(int i=0;i<v.size();i++){
		length=0;
		for(int j=0;j<v[0].size();j++)
			length+=pow(v[i][j],2);
		length=sqrt(length);
		for(int j=0;j<v[0].size();j++)
			v[i][j]/=length;
	}
	v=trans(v);
}

void orthonormalize_col(two_dvector &v){
	v=trans(v);
	for(int i=1;i<v.size();i++)
		for(int j=i-1;j>=0;j--)
			v[i]=v[i]-((v[j]*v[i])/(v[j]*v[j]))*v[j];
	v=trans(v);
	normalize_col(v);
}



two_dvector inversev(const two_dvector &v){
	matrix temp(v.size(),v[0].size());
	temp.co=v;
	return temp.inverse().co;
}

bool matrix::is_square(){
	if(row==col)
		return true;
	return false;
}
bool matrix::is_symmetric(){
	if(this->co==this->transpose().co)
		return true;
	return false;
}

//////////////////*functions for matrix*////////////////////

int matrix::get_row(){
	return row;
}
int matrix::get_col(){
	return col;
}
///basic operation////

void matrix::lu(){
	up=co;
	low=resize(row,row);
	permu=resize(row,row);
	set_to_identity(low);
	set_to_identity(permu);        
	for(int p=0;p<row-1;p++){
	    for(int k=p+1;k<row;k++){
	        for(int count=1;up[p][p]==0&&count<row-p;count++){
	        	up[p].swap(up[p+count]);
	        	permu[p].swap(permu[p+count]);
	        	per++;
			}
			if(up[p][p]!=0){
				low[k][p]=up[k][p]/up[p][p];
	        	up[k]=up[k]-(low[k][p]*up[p]);
			}
			else
				low[k][p]=0;
		}
	}
	two_dvector PA=permu*co;
	up=PA;
	set_to_identity(low);
	for(int p=0;p<row;p++){
	    for(int k=p+1;k<row;k++){
	    	if(up[p][p]!=0){
	    		low[k][p]=up[k][p]/up[p][p];
	        	up[k]=up[k]-(low[k][p]*up[p]);
			}
			else
				low[k][p]=0;
		}
	}
}

void matrix::ldu(){
	d=resize(row,row);
	lu();
	for(int i=0;i<row;i++){
		d[i][i]=up[i][i];
		if(up[i][i]!=0)
			up[i]=up[i]/up[i][i];
	}
}
matrix matrix::transpose() const{
	matrix ans(col,row);
	for(int i=0;i<row;i++)
		for(int j=0;j<col;j++)
			ans.co[j][i]=co[i][j];
	return ans;	
}

matrix matrix::inverse(){
		if(!is_square()){
			cout<<"Invalid input(should be a square matrix)"<<endl;
			matrix ans(0,0);
			return ans;
		}
		if(det()==0){
			cout<<"Invalid input(not full rank)"<<endl;
			matrix ans(0,0);
			return ans;
		}
		
		double fac;
		two_dvector aug(row,vector<double>(2*row));
		for(int i=0;i<row;i++)
			for(int j=0;j<row;j++)
				aug[i][j]=co[i][j];
		for(int i=0;i<row;i++)
			for(int j=row;j<2*row;j++)
					aug[i][j]=(j-row==i);
		for(int p=0;p<row;p++){
			for(int k=p+1;k<row;k++){
				for(int c=1;aug[p][p]==0;c++)
					aug[p].swap(aug[p+c]);
					fac=aug[k][p]/aug[p][p];
					for(int i=0;i<2*row;i++)
						aug[k][i]=aug[k][i]-fac*aug[p][i];	
				}
		}
		for(int p=row-1;p>=0;p--){
			for(int k=p-1;k>=0;k--){
				fac=aug[k][p]/aug[p][p];
				for(int i=0;i<2*row;i++)
					aug[k][i]=aug[k][i]-fac*aug[p][i];
			}
		}
		for(int p=0;p<row;p++){
			fac=aug[p][p];
			for(int i=0;i<2*row;i++)
				aug[p][i]/=fac;
		}
		matrix ans(row,row);
		for(int i=0;i<row;i++)
			for(int j=row;j<2*row;j++)
				ans.co[i][j-row]=aug[i][j];
		return ans;
	/*
	inver=resize(row,col);
	two_dvector u_inverse=resize(row,col);
	set_to_identity(u_inverse);
	two_dvector d_inverse=d;
	two_dvector l_inverse=resize(row,row);
	set_to_identity(l_inverse);
	for(int i=0;i<up.size();i++)
		for(int j=0;j<up[0].size();j++)
			u_inverse[i][j]=-up[i][j];
	for(int i=0;i<d.size();i++)
		if(d[i][i]!=0)
			d_inverse[i][i]=1/d[i][i];		
	for(int i=0;i<low.size();i++)
		for(int j=0;j<low[0].size();j++)
			l_inverse[i][j]=-low[i][j];
	inver=u_inverse*d_inverse*l_inverse*permu;
	return matrix(inver);
	*/
}

matrix matrix::get_rref(){
	double fac;
	matrix ans(row,col);
	rref=co;
	/*rref is U*/
	for(int p=0;p<row-1;p++){
	    for(int k=p+1;fabs(rref[p][p])!=0&&k<row;k++){
	        fac=rref[k][p]/rref[p][p];
	        rref[k]=rref[k]-fac*rref[p];
		}
	}
	/*remove the small number*/
	for(auto& r:rref)
		for(auto& c:r)
			if(fabs(c)<1e-6)
				c=0;
	/*number under pivots are zero*/
	for(int p=0;p<row;p++)
		for(int i=0;i<col;i++){
			if(fabs(rref[p][i])!=0){
				for(int k=p+1;k<row;k++){
					fac=rref[k][i]/rref[p][i];
					rref[k]=rref[k]-fac*rref[p];
				}
				break;
			}
	}
	for(int p=row-1;p>0;p--){
		for(int i=0;i<col;i++){
			if(rref[p][i]!=0){
				for(int k=p-1;k>=0;k--){
					fac=rref[k][i]/rref[p][i];
					rref[k]=rref[k]-fac*rref[p];
				}
			break;	
			}
		}
	}
	/*normalize*/
	for(int i=0;i<row;i++){		
		for(int j=0;j<col;j++){
			if(fabs(rref[i][j])!=0){
				for(int k=j+1;k<col;k++)
					rref[i][k]/=rref[i][j];
				rref[i][j]=1;
				break;
			}
		}
	}
	/*permutation*/
	for(int r=0;r<row;r++)
		for(int c=r+1;c<row;c++){
			if(fabs(rref[r][r])==0&&fabs(rref[c][r])!=0)
				rref[c].swap(rref[r]);
		}
	ans.co=rref;
	return ans;
}

int matrix::rank(){
	get_rref();
	r=0;
	for(int i=0;i<row;i++){
		for(int j=i;j<col;j++){
			if(fabs(rref[i][j])!=0){
				r++;
				break;
			}	
		}
	}
	return r;
}

matrix matrix::project(){
	if(rank()==col){
		matrix ans=(*this)*((transpose()*(*this)).inverse())*transpose();
		return ans;
	}
	else{
		cout<<"Invalid Input(not full column rank)"<<endl;
		matrix ans(0,0);
		return ans;
	}
}

double matrix::det(){
	if(!is_square()){
		cout<<"Invalid input(should be a square matrix)"<<endl;
		return 0;
	}	
	ldu();
	determine=1;
	for(int i=0;i<row;i++)
		determine*=d[i][i];
	if(per%2==1)
		determine=-determine;
	return determine;
}


/*four subspaces and orthogonaliztion*/

two_dvector matrix::rowspace(){
	row_space.clear();
	for(int i=0;i<rank();i++)
		row_space.push_back(rref[i]);
	return row_space;
}
two_dvector matrix::colspace(){
	pivot_col.clear();
	col_space.clear();
	rank();
	for(int i=0;i<row;i++){
		for(int j=i;j<col;j++){
			if(rref[i][j]!=0){
				pivot_col.push_back(j);
				break;
			}
		}
	}
	for(int i=0;i<r;i++)
		col_space.push_back(transpose().co[pivot_col[i]]);
	col_space=trans(col_space);
	return col_space;
}
two_dvector matrix::nullspace(){
	pivot_col.clear();
	free_col.clear();
	null_space.clear();
	if(col==rank()){
		return null_space;
	}	
	for(int i=0;i<row;i++){
		for(int j=i;j<col;j++){
			if(fabs(rref[i][j])!=0){
				pivot_col.push_back(j);
				break;
			}
		}
	}
	null_space=resize(col,col-r);
	two_dvector tempr=rref;
	for(int i=row-1;i>=0;i--){
		for(int j=r-1;j>=0;j--){
			tempr[i].erase(tempr[i].begin()+pivot_col[j]);
		}
	}
	for(int i=0;i<r;i++)
		for(int j=0;j<col-r;j++)
			null_space[i][j]=-tempr[i][j];
	for(int i=r;i<col;i++)
		for(int j=0;j<col-r;j++)
			null_space[i][j]=(j==i-r);
	for(int i=0;i<col;i++)
		if(!is_in(pivot_col,i))
			free_col.push_back(i);
	tempr=null_space;
	
	for(int i=0;i<r;i++)
		null_space[pivot_col[i]]=tempr[i];
	for(int i=0;i<col-r;i++)
		null_space[free_col[i]]=tempr[i+r];
	return null_space;
}
two_dvector matrix::leftnullspace(){
	left_null_space=transpose().nullspace();
	return trans(left_null_space);
}
void matrix::GQR(){
	long double cos=0,sin=0;
	GQ.resize(row);
	for(int i=0;i<GQ.size();i++)
		GQ[i].resize(row,0);
	set_to_identity(GQ);
	two_dvector tempQ=GQ;
	GR=co;
	for(int i=0;i<col-1;i++){
		for(int j=i+1;j<row;j++){
			givens_rotation(GR[i][i],-GR[j][i],cos,sin);
			set_to_identity(tempQ);
			tempQ[i][i]=cos;
			tempQ[j][j]=cos;
			tempQ[j][i]=sin;
			tempQ[i][j]=-sin;
			GR=tempQ*GR;
			GQ=GQ*trans(tempQ);
		}
	}
}
void matrix::GSQR(){
	GSQ=trans(co);
	GSR=resize(col,col);
	for(int i=1;i<col;i++)
		for(int j=i-1;j>=0;j--)
			GSQ[i]=GSQ[i]-((GSQ[j]*GSQ[i])/(GSQ[j]*GSQ[j]))*GSQ[j];
	GSQ=trans(GSQ);
	normalize_col(GSQ);
	for(int i=0;i<col;i++)
		for(int j=0;j<col;j++)
			GSR[i][j]=trans(GSQ)[i]*trans(co)[j];
}

/*eigenvalue and diagonalization*/

set< double , greater<double> > matrix::eigenvalue(){
	evalue.clear();
	evaluev.clear();
	if(!is_square()){
		cout<<"Invalid input(should be a square matrix)"<<endl;
		return evalue;
	}
	if(!is_symmetric())
		cout<<"Warning!!!This matrix isn't symmetric, may have complex eigenvalues and eigenvectors.\n\n";
	matrix old_dia=*this;
	matrix new_dia=*this;
	two_dvector tempQ=co;
	two_dvector tempR=co;
	Q=co;
	set_to_identity(Q);
	bool thres=true;
	while(thres){
		old_dia.GQR();
		new_dia.co=(old_dia.GR)*(old_dia.GQ);
		Q=Q*old_dia.GQ;
		thres=false;
		for(int i=0;i<row;i++)
			if(fabs(new_dia.co[i][i]-old_dia.co[i][i])>1e-6)
				thres=true;
		old_dia=new_dia;
	}
	for(int i=0;i<row;i++)
		evalue.insert(round(1e7*new_dia.co[i][i])/1e7);
	for(int i=0;i<row;i++)
		evaluev.push_back(round(1e7*new_dia.co[i][i])/1e7);
	sort_for_svd();	
	return evalue;
}
two_dvector matrix::eigenvector(){
	evector.clear();
	eigenvalue();
	if(!is_square()){
		return evector;
	}
	if(is_symmetric()){
		Q_lambda_QT();
		/*streching so that fabs>=1(better to read)*/
		for(int i=0;i<col;i++){
			double smallest=fabs(Q[0][i]);
			for(int j=1;j<row;j++){
				if(fabs(Q[j][i])<smallest)
					smallest=fabs(Q[j][i]);
			}
			for(int j=0;j<row;j++){
				Q[j][i]*=1./smallest;
			}
		}
		/*round off some number*/
		for(auto& r:Q)
			for(auto& c:r)
				c=round(c*1e3)/1e3;
		evector=Q;
		return evector;
	}
	/*if not symmetric*/
	/*have a matrix type of Identity*/
	matrix I(row,col);
	set_to_identity(I.co);
	/*compute the nullspace(dangerous here)*/
	for(auto eval:evalue){
		matrix tempm=*this-eval*(I);
		two_dvector tempv=tempm.nullspace();
		for(auto null:trans(tempv)){
			if(trans(tempv).size()!=0)/*nullity*/
				evector.push_back(null);
			else{
				cout<<"Round off error!!!"<<endl;
				evector.clear();
				return evector;
			}
		}
	}
	evector=trans(evector);
	return evector;
}

void matrix::Q_lambda_QT(){
	Q.clear();
	lambda.clear();
	eigenvalue();
	lambda=resize(row,row);
	for(int i=0;i<row;i++)
		lambda[i][i]=evaluev[i];
}

void matrix::S_lambda_S1(){
	two_dvector ev=trans(evector);
	for(int i=0;i<ev.size();i++)
		S.push_back(ev[i]);
	for(int i=0;i<row;i++)
		lambda[i][i]=evaluev[i];
}

void matrix::svd(){
	matrix temp=transpose()*(*this);
	temp.Q_lambda_QT();
	sigma=resize(row,col);
	for(int i=0;i<rank();i++)
		sigma[i][i]=sqrt(temp.lambda[i][i]);
	V=temp.Q;
	/*
	temp=(*this)*transpose();
	temp.Q_lambda_QT();
	U=temp.Q;
	*/
	U=resize(row,row);
	for(int i=0;i<r;i++)
		U[i]=(1/sigma[i][i])*(co*trans(V)[i]);
	for(int i=r;i<row;i++)
		U[i]=leftnullspace()[i-r];
	U=trans(U);
	orthonormalize_col(U);
}

/*Pseudo inverse*/
matrix matrix::pseudo_inverse(){
	
	if(row==col&&row==rank()){
		return inverse();
	}
	
	svd();
	matrix ans(col,row);
	two_dvector sigma_plus=resize(col,row);
	for(int i=0;i<min(row,col);i++)
		if(sigma[i][i]!=0)
			sigma_plus[i][i]=1/sigma[i][i];
	ans.co=V*sigma_plus*trans(U);
	return ans;
	
} 
/*approximation*/

matrix matrix::rank_n_appr(int n){
	if(n<rank()){
		matrix ans(row,col);
		two_dvector temp=resize(row,col);
		svd();
		for(int i=0;i<n;i++){
			for(int j=0;j<row;j++)
				temp[j]=U[j][i]*trans(V)[i];
			ans.co=ans.co+sigma[i][i]*temp;
		}
		return ans;
	}	
	else{
		cout<<"Invalid input, n should < "<<r<<endl;
		return matrix(0,0);
	}	
}
/*Ax=b*/

matrix Axb(matrix& a,matrix& b){
	if(a.row==b.row){
		if(a.rank()==a.row)
			return a.inverse()*b; 
		return a.pseudo_inverse()*b;
	}
	else{
		cout<<"Invalid input"<<endl;
		return matrix(0,0);
	}
}

/* condition number */
double matrix::con_num(){
	if(rank()<min(row,col)){
		return 1e20;
	}
	else{
		svd();
		return (sigma[0][0]/sigma[row-1][col-1]);
	}
}
/*helper function*/

void matrix::sort_for_svd(){
	Q=trans(Q);
	for(int i=0;i<evaluev.size()-1;i++)/*bubble sort*/
		for(int j=0;j<evaluev.size()-i-1;j++)
			if(evaluev[j]<evaluev[j+1]){
				swap(evaluev[j],evaluev[j+1]);
				Q[j].swap(Q[j+1]);
			}
	Q=trans(Q);	
}

/////////////////*get functions*//////////////////

void matrix::getlu(ostream& outs){
	lu();	
	if(per)
		outs<<"P= "<<endl<<permu<<endl;
	outs<<"A= "<<endl<<co<<endl;
	outs<<"L= "<<endl<<low<<endl;
	outs<<"U= "<<endl<<up<<endl;
}

void matrix::getldu(ostream& outs){
	ldu();
	if(per)
		outs<<"P= "<<endl<<permu<<endl;
	outs<<"A= "<<endl<<co<<endl;
	outs<<"L= "<<endl<<low<<endl;
	outs<<"D= "<<endl<<d<<endl;
	outs<<"U= "<<endl<<up<<endl;
}
void matrix::getGQR(ostream& outs){
	GQR();
	cout<<"Q= "<<endl<<GQ<<endl;
	cout<<"R= "<<endl<<GR<<endl;
}

void matrix::getGSQR(ostream& outs){
	GSQR();
	cout<<"Q= "<<endl<<GSQ<<endl;
	cout<<"R= "<<endl<<GSR<<endl;
}

void matrix::getQ_lambda_QT(ostream& outs){
	if(!is_square()){
		outs<<"Invalid input(Not square)"<<endl;
		return;
	}
	if(!is_symmetric()){
		outs<<"Invalid input(Not symmetric)"<<endl;
		return ;
	}
	Q_lambda_QT();	
	outs<<"Q="<<endl<<Q<<endl;
	outs<<"¢GN="<<endl<<lambda<<endl;
	outs<<"QT="<<endl<<trans(Q)<<endl;
}

void matrix::getS_lambda_S1(ostream& outs){
	if(!is_square())
		return;
	S.clear();
	lambda.clear();
	eigenvector();
	lambda=resize(row,row);
	if(evector[0].size()!=row){
		outs<<"Not diagonalizable(Not enough eigenvector)"<<endl;
		return;
	}
	S_lambda_S1();		
	outs<<"S="<<endl<<S<<endl;
	outs<<"¢GN="<<endl<<lambda<<endl;
	outs<<"S-1="<<endl<<inversev(S)<<endl;
}

void matrix::getsvd(ostream& outs){
	svd();
	outs<<"U="<<endl<<U<<endl;
	outs<<"¢GU="<<endl<<sigma<<endl;
	outs<<"VT="<<endl<<trans(V)<<endl;
}


bool is_in(const vector<int>& v,int ele){
	for(auto x :v)
		if(x==ele)
			return true;
	return false;
}

int sgn(long double x){
	if(x>=0)
		return 1;
	return -1;
}

void givens_rotation(double a,double b,long double &c,long double &s){
	if(b==0){
		c=sgn(a);
        s=0;
	}
    else if(a==0){
    	c=0;
        s=sgn(b);
	}  
    else if(fabs(a)>fabs(b)){
    	long double tan=b/a;
        long double u=sgn(a)*sqrt(1+tan*tan);
        c=1/u;
        s=c*tan;
	}  
    else{
    	long double tan=a/b;
        long double u=sgn(b)*sqrt(1+tan*tan);
        s=1/u;
        c=s*tan;
	}
}

/*destructor*/
matrix::~matrix(){
}

void print_manu(ostream& outs){
	outs<<"\n1:addition / 2:subtraction / 3:multication / 4:transpose  / 5:LU decomposition\n\n";
	outs<<"6:LDU decomposition / 7:inverse / 8:rank / 9:row reduced echelon form / 10:Column space\n\n";
	outs<<"11:Null space / 12:Row space / 13:Left null space / 14:QR decomposition(Givens rotation)\n\n";
	outs<<"15:QR decomposition(Gram-Schmidt method) / 16:projection / 17:determinant / 18:Eigenvalues\n\n";
	outs<<"19:Eigenvectors(If not symmetric due to numerical instability, this may have errors)\n\n";
	outs<<"20:Q ¢GN QT diagonalization(for symmetric matrix)\n\n";
	outs<<"21:S ¢GN S-1 diagonalization(Due to the numerical instability , this may have errors)\n\n";
	outs<<"22:Singular value decomposition / 23:Best rank n approximation / 24:Pseudo inverse\n\n";
	outs<<"25:Solving Ax=b 26:Condition number\n\n";
	outs<<"Press 0 to quit"<<endl;
}

int main(){
	cout<<"-----------Welcome-----------"<<endl;
	int select,n;
	print_manu(cout);
	matrix a,b;
	while(cin>>select){
		if(select==0)
			break;
		switch(select){
			case 1:
				cin>>a>>b;
				cout<<a+b;
				break;
			case 2:
				cin>>a>>b;
				cout<<a-b;
				break;
			case 3:
				cin>>a>>b;
				cout<<a*b;
				break;
			case 4:
				cin>>a;
				cout<<a.transpose();
				break;
			case 5:
				cin>>a;
				a.getlu(cout);
				break;
			case 6:
				cin>>a;
				a.getldu(cout);
				break;
			case 7:
				cin>>a;
				cout<<a.inverse();
				break;
			case 8:
				cin>>a;
				cout<<"rank="<<a.rank()<<endl;
				break;
			case 9:
				cin>>a;
				cout<<"rref="<<endl<<a.get_rref()<<endl;
				break;
			case 10:
				cin>>a;
				cout<<"The column space is"<<endl<<"span"<<endl<<a.colspace()<<endl;
				break;
			case 11:
				cin>>a;
				cout<<"The nullspace is ";
				if(a.rank()==a.get_col())
					cout<<"Z"<<endl;
				else
					cout<<endl<<"span"<<endl<<a.nullspace()<<endl;
				break;
			case 12:
				cin>>a;
				cout<<"The row space is "<<endl<<"span"<<endl<<a.rowspace()<<endl;
				break;
			case 13:
				cin>>a;
				cout<<"The left nullspace is ";
				if(a.get_row()==a.rank())
					cout<<"Z"<<endl;
				else
					cout<<endl<<"span"<<endl<<a.leftnullspace()<<endl;
				break;
			case 14:
				cin>>a;
				a.getGQR(cout);
				break;
			case 15:
				cin>>a;
				a.getGSQR(cout);
				break;
			case 16:
				cin>>a;
				cout<<a.project();
				break;
			case 17:
				cin>>a;
				cout<<"det="<<a.det()<<endl;
				break;
			case 18:
				cin>>a;
				cout<<"The eigenvalues are "<<a.eigenvalue()<<endl;
				break;
			case 19:
				cin>>a;
				cout<<"The eigenvectors are "<<endl<<a.eigenvector()<<endl;
				break;
			case 20:
				cin>>a;
				a.getQ_lambda_QT(cout);
				break;
			case 21:
				cin>>a;
				a.getS_lambda_S1(cout);
				break;
			case 22:
				cin>>a;
				a.getsvd(cout);
				break;
			case 23:
				cin>>a;
				cout<<"Please input the rank you want to approximate"<<endl;
				cin>>n;
				cout<<a.rank_n_appr(n);
				break;
			case 24:
				cin>>a;
				cout<<a.pseudo_inverse()<<endl;
				break;
			case 25:
				cin>>a>>b;
				cout<<"x="<<endl<<Axb(a,b);
				break;
			case 26:
				cin>>a;
				cout<<"condition number= "<<a.con_num()<<endl;
				break; 
		}
		print_manu(cout);
	}
	cout<<"----------Thank you-------------"<<endl;
	return 0;
}

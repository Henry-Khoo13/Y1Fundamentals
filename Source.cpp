#include <windows.h>
#include <tchar.h>

#include <iostream>
#include <string>

#include<math.h> //For the Pow (Power)

using namespace std;
namespace MatrixThings {
	struct myMat {			// allows for matrices up to size 4*4
		int numRows;		// number of rows
		int numCols;
		int data[16];		// data are stored in row order
	};

	myMat zeroMat(int r, int c) {
		// create a matrix with r rows and c columns, filled with zeros
		myMat m;			// define matrix
		m.numRows = r;		// set size
		m.numCols = c;
		for (int ct = 0; ct < 16; ct++) m.data[ct] = 0;	// set elems to 0
		return m;			// return matrix
	}

	int getElem(myMat m, int r, int c) {
		// return the item at m[r,c]   where r is 0..m.numRows-1 etc
		return m.data[r * m.numCols + c];
	}

	void setElem(myMat& m, int r, int c, int val) {
		// set element m[r,c] to val
		m.data[r * m.numCols + c] = val;
	}

	myMat mFromStr(string s) {
		// create a matrix from string s
		// string of form "1,2;3,4;5,6"   ; separates rows and , separates columns ... No error check
		int ms;
		if (s.length() > 0) ms = 1; else ms = 0;
		myMat m = zeroMat(ms, ms);						// is s empty create 0*0 matrix, else start with 1*1 matrix
		int mndx = 0;									// used to index into array
		string sub = "";								// sub string - numbers between , or ; set empty
		for (int ct = 0; ct < s.length(); ct++) {		// each char in turn
			if ((s[ct] == ',') || (s[ct] == ';')) {	// if , or ; then sub contains number
				m.data[mndx++] = stoi(sub);				// convert number to integer, put in data array
				sub = "";								// clear sub string
				if (s[ct] == ';') m.numRows++;			// if found ; indicates an extra row
				else if (m.numRows == 1) m.numCols++;	// if , then (if in row 1) increase count of columns
			}
			else sub = sub + s[ct];						// add character to sub string
		}
		if (sub.length() > 0) m.data[mndx++] = stoi(sub);// add last sub string
		return m;
	}

	myMat mGetRow(myMat m, int row) {
		// create a matrix from m, having one row
		myMat res = zeroMat(1, m.numCols);		// create a matrix of 1 row
		for (int col = 0; col < m.numCols; col++)		// for each element in row
			res.data[col] = getElem(m, row, col);		// copy col element to res
		return res;
	}

	myMat mGetCol(myMat m, int col) {
		// create a matrix from m, having one col
		myMat res = zeroMat(m.numRows, 1);		// create a matrix of 1 column
		for (int row = 0; row < m.numRows; row++)		// for each element in col
			res.data[row] = getElem(m, row, col);		// copy row element to res
		return res; //Return Matrix
	}
	
	int dotProd(myMat v1, myMat v2) { //This function takes two matrices.
		int f = 0;//F is the total of the dot product and will be the result returned.
		for (int a = 0; a < v1.numRows; a++) //The code processes every row of both
		{
			for (int b = 0; b < v2.numCols; b++) //And processes Every Column in every row a.k.a getting every element of each row
			{
				int c = getElem(v1, a, b);//Using the get element function to return the value in matrix 1
				int d = getElem(v2, a, b);//Using the get element function to return the value in matrix 2 
				//It gets the same element of both matrices that the program is current on. (The column and row location is the same)
				int e = c * d;//Multiply both gathered elements together
				f = f + e;//Add the result to a total
			}
		}
		return f;//Return the total/result
	}
	
	myMat mTranspose(myMat m) {//This function takes one matrix.
		myMat res = zeroMat(m.numCols, m.numRows);//Inversed Number of Columns and Rows 3*2 -> 2*3
		//Now to go throuh each element of the current matrix.
		for (int a = 0; a < m.numRows; a++)//Processing Every Row..
		{
			for (int b = 0; b < m.numCols; b++)//Processing Every Column in every row a.k.a getting every element of each row
			{
				int c = getElem(m, a, b);//Getting the Current Element
				setElem(res, b, a, c);//Putting it in the inversed matrix [1,2] -> [2,1] (Reversed Column and Row)
			}
		}
		return res;//Output the new matrix.
	}
	
	myMat mAdd(myMat m1, myMat m2) {//Assume both arrays are the same size in both dimensions
		myMat res = zeroMat(m1.numRows, m2.numCols);//Creates an empty matrix using the sizes of the two other matrices.
		//It doesn't matter which rows or cols are used since they should be both the same size.
		for (int a = 0; a < m1.numRows; a++)//Processing Every Row..
		{
			for (int b = 0; b < m2.numCols; b++)//Processing Every Column in every row a.k.a getting every element of each row
			{
				int c = getElem(m1, a, b);//Using the get element function to return the value in matrix 1
				int d = getElem(m2, a, b);//Using the get element function to return the value in matrix 2 
				int Total = c + d;//Using the two gathered elements, they are totalled together
				setElem(res, a, b, Total);//The total of the two elements is inserted into the same column and row location in the empty matrix.
			}
		}
		return res;//The now filled "Empty" matrix is returned.
	}
	
	myMat mMult(myMat m1, myMat m2) {//Takes two matrices
		myMat FullRes = zeroMat(m1.numRows, m2.numCols);//Uses the number of rows from the first matrix and number of columns from the second matrix.
		for (int m1Row = 0; m1Row < m1.numRows; m1Row++)//Goes through each element of the Full res (Each Row)
		{
			for (int m2Column = 0; m2Column < m2.numCols; m2Column++)//Goes through each element of the Full res  (And Each Column to get each element)
			{
				myMat CurrentRow1 = mGetRow(m1, m1Row);//Get the Row of the current Element
				myMat CurrentCol2 = mGetCol(m2, m2Column);//Get the Column of the current Element
				int Total = 0;
				for (int CR1Cols = 0; CR1Cols < CurrentRow1.numCols; CR1Cols++) {//Goes through each element of both Column and Row matrix 
					int c = getElem(CurrentRow1, 0, CR1Cols);//Gets Element from Row
					int d = getElem(CurrentCol2, CR1Cols, 0);//Gets Element from column
					int e = c * d;//Multiples 
					Total = Total + e;//Sums up the total of all the elemts in both the Row and Column Matrix
					//cout << " C " << c << " D " << d << " T " << Total << endl;
				}
				setElem(FullRes, m1Row, m2Column, Total);//Sets the result into the result matrix.
			}
		}
		return FullRes;//Returns the result. 
	}

	myMat mScalar(myMat m, int i) {//Takes a Matrix and the Scalar that will be used on the matrix.
		myMat Res = zeroMat(m.numRows, m.numCols);//Creates a Matrix of the same size as the original.
		for (int mRow = 0; mRow < m.numRows; mRow++)//Goes through each row of the original
		{
			for (int mColumn = 0; mColumn < m.numCols; mColumn++)//Goes through each column of the original (Thus element)
			{
				int c = getElem(m, mRow, mColumn);//Gets the element from the original
				setElem(Res, mRow, mColumn, c * i);//Multiplies it againt the scalar then places it into the new result matrix.
			}
		}
		return Res;//Returns the result matrix.
	}

	myMat mSubM(myMat m, int i, int j) {//Small Matrix elinating row and column in parameters (Used for determinant)
		myMat Res = zeroMat(m.numRows - 1, m.numCols - 1);//A matrix is created which has one less row and column than than the original
		int RowPos = 0, ColPos = 0, mRow, mColumn;//The Loop and Current position that the result matrix is being filled are here (Loop ones needed to compare)
		for (mRow = 0; mRow < m.numRows; mRow++)//Each Row
		{
			for (mColumn = 0; mColumn < m.numCols; mColumn++)//And Each column of the original
			{
				if (mColumn == j ) {}//We check if it's the column we want to eleminate, if so the code does nothing
				else if (mRow == i) {}//We check if it's the row we want to eleminate, if so the code does nothing
				else { //Otherwise...
					int c = getElem(m, mRow, mColumn);//The code gets the current element 
					setElem(Res, RowPos, ColPos, c);//Places the element into the result matrix using the current pos variables made before.
					ColPos = ColPos + 1;//Current pos of the column is updated.
				}
			}
			ColPos = 0;
			if (mColumn == j) {}//We check if it's the column we want to eleminate, if so the code does nothing
			else if (mRow == i) {}//We check if it's the row we want to eleminate, if so the code does nothing
			else RowPos = RowPos + 1;//Current pos of the Row is updated.
		}
		return Res;//The result matrix is then outputted.
	}

	int mDeterminant(myMat m) {//Takes a Matrix
		int i = 1;
		if (m.numRows == 1) {//Checks if it only has one Row 1*1 matrix
			return getElem(m, 0, 0);//Outputs the first element otherwise
		}
		else if (m.numRows == 2) {//Checks if has two rows 2*2 matrix
			return (getElem(m, 0, 0) * getElem(m, 1, 1)) - (getElem(m, 0, 1) * getElem(m, 1, 0)); //Returns the determinant after some cals
		}
		else {
			int ans = 0;//Answer summing each cofactor 
			for (int ct = 0; ct < m.numCols; ct++) ans += pow(-1, ct + 2) * getElem(m, 0, ct) * mDeterminant(mSubM(m, 0, ct));//Uses fibonacci to find the determinant
			return ans;//returns answer
		}
		return i;//Returns 1 if error.
	}

	myMat mAdjoint(myMat m) {
		myMat Res = zeroMat(m.numRows, m.numCols);
		for (int mRow = 0; mRow < m.numRows; mRow++)
		{
			for (int mColumn = 0; mColumn < m.numCols; mColumn++)
			{
				int C = pow(-1, mRow + mColumn) * mDeterminant(mSubM(m, mRow, mColumn));
				setElem(Res, mColumn, mRow, C);
			}
		}
		return Res;
	}

	myMat mCramer(myMat A, myMat b) {//Takes two matrices.
		int d = mDeterminant(A);//Get the determinant with previous determinant function
		myMat m, x;//Creates an M and X matrix 
		x = zeroMat(A.numRows, A.numCols);//Initalises X matrix with empty amount of rows and cols
		for (int c = 0; c < b.numCols+1;c++) {//Goes through each col
			m = A;//Sets Matrix m to A
			for (int r = 0; r < b.numRows+1; r++) {// Goes through each row
				int BR = getElem(b, r, 0);//Gets each element from b 
				setElem(m, r, c, BR);//Sets the result into M
				
			}
			int e = mDeterminant(m) / d;//Divides the determinant of M by the determinant of A
			setElem(x, c, 0, e);//Sets the answer into the result matrix X
		}
		return x;//Returns the result matrix X.
	}


//Scrapped Code
/*
	myMat mInverse(myMat m) {
		myMat Res = zeroMat(m.numRows, m.numCols);
		int det = mDeterminant(m);
		//cout << "Det: " << det;
		float DivDet = 1 / det;
		for (int mRow = 0; mRow < m.numRows; mRow++)
		{
			for (int mColumn = 0; mColumn < m.numCols; mColumn++)
			{
				int c = getElem(m, mRow, mColumn);
				int d = (1 / det) * c;
				cout << "Value C : " << c << endl;
				cout << "Value D : " << d << endl;
				setElem(Res, mRow, mColumn, d);
			}
		}
		return Res;
	}
	*/
/*
	for (int RowV1Length = 0; RowV1Length < res.numCols; RowV1Length++)
	{
		for (int ColumnV2Length = 0; ColumnV2Length < res2.numCols; ColumnV2Length++)
		{
			int Val1 = res.data[RowV1Length];
			int Val2 = res2.data[ColumnV2Length];
			int MultiplyRes = Val1 * Val2;
			cout << "Val1: " << Val1 << endl;
			cout << "Val2: " << Val2 << endl;
			cout << "MultiplyRes: " << MultiplyRes << endl;
			FullRes.data[MultiplyRes];//Not sure this works.

		}//Currently getting each value from row 1 of Matrix A and column 1 of Matrix B and multiplying Will place into new matrix
	}
	*/
/*Scrapped Code
					myMat m1RowRes = mGetRow(m1, m1Row);
				myMat m2ColRes = mGetCol(m2, m2Column);
				for (int m1RowResLength = 0; m1RowResLength < m1RowRes.numCols; m1RowResLength++) {
					int Total = 0;
					for (int m2ColResLength = 0; m2ColResLength < m2ColRes.numRows; m2ColResLength++) {
						int MultiRes = m1RowRes.data[m1RowResLength] * m2ColRes.data[m2ColResLength];
						cout << "Row: " << m1Row << endl << "Row Num: " << m1RowResLength << endl << "Row Elm:" << m1RowRes.data[m1RowResLength] << endl;
						cout << "Col: " << m2Column << endl << "Col Num: " << m2ColResLength << endl << "Col Elm:" << m2ColRes.data[m2ColResLength] << endl;
						cout << "Multi Res: " << MultiRes << endl;
						Total = Total + MultiRes;
						cout << "Total Res: " << Total << endl;
						cout << "\n";
						setElem(FullRes, m1Row, m2Column, Total);
					}
				}
	
	*/
/*
for (int e = 0; e < res.numRows; e++) {//Processing Every Row
	for (int f = 0; f < res.numCols; f++)//Processing Every Column in every row a.k.a getting every element of each row
	{
		int g = getElem(v1, e, f);//Using the get element function to return the value in matrix res
		Total2 += g;//Adding each element to a new totalling variable to collect the sum of all values in matrix res.
	}
}
*/
/*				int c = getElem(m, mRow, mColumn);
			//cout << "C : " << c;
			//cout << "R: " << RowPos << " C: " << ColPos << endl;
			cout << "R: " << RowPos << " C: " << ColPos << endl;
			cout << "Test";
			cout << "R: " << mRow << " C: " << mColumn << endl;
			ColPos = ColPos + 1;
			if (mRow == i) {
				cout << "R: " << RowPos << " C: " << ColPos << endl;
				cout << "";
				//cout << "R: " << mRow << " C: " << mColumn << endl;
			}
			else {
				//setElem(Res, RowPos, ColPos, c);
				//ColPos = ColPos + 1;
			}
			//ColPos = ColPos + 1;
						RowPos = RowPos + 1;
		if (mRow == i) {}
		else {
			//RowPos = RowPos + 1;
		}
*/
/*		m8 = mFromStr("91");
		printMat("m8", m8);
		DetTest = mDeterminant(m8);
		cout << "Determinant: " << DetTest << endl;

		m9 = mFromStr("5,2; 3,4");
		printMat("m9", m9);
		DetTest = mDeterminant(m9);
		cout << "Determinant: " << DetTest << endl;

		m11 = mFromStr("2,2,0;-2,1,1;3,0,1");
		printMat("m10", m11);
		DetTest = mDeterminant(m11);
		cout << "Determinant: " << DetTest << endl;

		//Testing Shrunk Matrix
		m12 = mSubM(m11, 0, 0);
		printMat("m12", m12);
		DetTest = mDeterminant(m12);

		cout << "Cramer" << endl;
		m10 = mFromStr("16; 18");
		m13 = mCramer(m9,m10);
		printMat("m13", m13);

		//m12 = mInverse(m10);
		//printMat("m12", m12);

		m1 = mFromStr("1,2,3;4,-5,6");
		m2 = mFromStr("2,-5;3,7;1,4");
		printMat("m1", m1);
		printMat("m2", m2);
		//Transpose
		cout << "Transpose" << endl;
		m3 = mTranspose(m1);
		printMat("m3", m3);
		//Multi
		cout << "Multiply" << endl;
		m4 = mMult(m2, m1);
		printMat("m4", m4);
		//Dot Product
		cout << "Dot Product" << endl;
		//m6 = mFromStr("1,2;-5,6");
		//m7 = mFromStr("2,-5;3,7");
		m6 = mFromStr("1;1");
		m7 = mFromStr("-1;0");
		printMat("m6", m6);
		printMat("m7", m7);
		DotProductRes = dotProd(m6, m7);
		cout << "Dot Product: " << DotProductRes << endl;
		//https://www.theclickreader.com/dot-products-and-matrix-multiplication/
	
		//Add is functional
		m3 = mFromStr("1,2;-5,6");
		m4 = mFromStr("2,-5;3,7");
		printMat("m3", m3);
		printMat("m4", m4);
		m5 = mAdd(m3, m4);
		printMat("m5", m5);
		//m3 = mGetCol(m1, 2);
		//printMat("Row1_m1", m3);
		//m5 = mMult(m2, m1);
		//printMat("m2*m1", m5);
		//myMat res = mGetRow(m1, 0);
		//printMat("m1.col.res.test", res);
		//myMat res2 = mGetCol(m1, 0);
		//printMat("m1.Col.res.test", res2);
		//int elm0 = res.data[0];
		//cout << "Element0 Output Test: " << elm0 << endl;
		//int elm1 = res.data[1];
		//cout << "Element1 Output Test: " << elm1 << endl;
		//int elm2 = res.data[2];
		//cout << "Element1 Output Test: " << elm2 << endl;
		//return 0;*/
/*
//https://www.mathsisfun.com/algebra/matrix-multiplying.html
//Go through the new matrix and total up all the elements within it and add them together.
*/

	void printMat(const char* mess, myMat m) {
		cout << mess << " = " << "\n";
		for (int r = 0; r < m.numRows; r++) {
			for (int c = 0; c < m.numCols; c++)
				cout << getElem(m, r, c) << "\t";
			cout << "\n";
		}
		cout << "\n";
	}
	
	int UsingMatrix() {
		cout << "Henry's Matrix Program\n";
		cout << "Elements Taken from: Ric's Matrix Example Program\n";
		myMat m1, m2, m3, m4, m5, m6, m7,m8,m9,m10,m11,m12,m13;
		
		m1 = mFromStr("1,2,3;4,-5,6");
		m2 = mFromStr("2,-5;3,7;1,4");

		printMat("m1", m1);
		printMat("m2", m2);
		cout << "\n" << endl;

		cout << "Outputting Column 0 of m1 " << endl;
		m3 = mGetCol(m1, 0);
		printMat("m3",m3);

		cout << "Outputting Column 1 of m1 " << endl;
		m3 = mGetCol(m1, 1);
		printMat("m3", m3);

		cout << "Outputting Column 2 of m1 " << endl;
		m3 = mGetCol(m1, 2);
		printMat("m3", m3);

		cout << "Dot Product of m1 and m1" << endl;
		int DotProductRes = dotProd(m1, m1);
		cout << "Dot Product Result: " << DotProductRes << endl;
		cout << "\n" << endl;

		cout << "Dot Product of m2 and m2" << endl;
		DotProductRes = dotProd(m2, m2);
		cout << "Dot Product Result: " << DotProductRes << endl;
		cout << "\n" << endl;

		cout << "Transpose m1" << endl;
		m3 = mTranspose(m1);
		printMat("m3", m3);
		cout << "\n" << endl;

		cout << "Transpose m2" << endl;
		m3 = mTranspose(m2);
		printMat("m3", m3);
		cout << "\n" << endl;

		cout << "Addition, M1 + M1" << endl;
		m3 = mAdd(m1, m1);
		printMat("m3", m3);
		cout << "\n" << endl;

		cout << "Addition, M2 + M2" << endl;
		m3 = mAdd(m2, m2);
		printMat("m3", m3);
		cout << "\n" << endl;

		cout << "Multiply m1 * m2" << endl;
		m3 = mMult(m1, m2);
		printMat("m3", m3);
		cout << "\n" << endl;


		cout << "Multiply m2 * m1" << endl;
		m3 = mMult(m2, m1);
		printMat("m3", m3);
		cout << "\n" << endl;

		cout << "Scalar m1 by 3" << endl;
		m3 = mScalar(m1, 3);
		printMat("m3", m3);
		cout << "\n" << endl;

		cout << "Scalar m2 by 2" << endl;
		m3 = mScalar(m2, 2);
		printMat("m3", m3);
		cout << "\n" << endl;

		cout << "m2 mSubM(0,0) Remove Column 0 and Row 0" << endl;
		m3 = mSubM(m2, 0,0);
		printMat("m3", m3);
		cout << "\n" << endl;

		m4 = mFromStr("91");
		printMat("m4", m4);
		int DetTest = mDeterminant(m4);
		cout << "Determinant of m4: " << DetTest << endl;

		m5 = mFromStr("5,2; 3,4");
		printMat("m5", m5);
		DetTest = mDeterminant(m5);
		cout << "Determinant of m5: " << DetTest << endl;

		m6 = mFromStr("2,2,0;-2,1,1;3,0,1");
		printMat("m6", m6);
		DetTest = mDeterminant(m6);
		cout << "Determinant of m6: " << DetTest << endl;

		cout << "Cramer : A = m7, b = m8 "<<endl;
		m7 = mFromStr("5,2;3,4");
		m8 = mFromStr("16;18");
		m9 = mCramer(m7, m8);
		printMat("m7", m7);
		printMat("m8", m8);
		printMat("m9", m9);

		return 0;
	}
}
namespace RecursionWork {
	int FibonacciRecursion(int n) {
		if (n <= 1) {
			return n; //Termination Condition (Outputs the Result)
		}
		else {
			return FibonacciRecursion(n - 1) + FibonacciRecursion(n - 2); //Fibonnacci Recurssion Code
		}
	}
	int FibonacciIterative(int n) {
		int fib, fib1, fib2;
		fib1 = 1;
		fib2 = 0;
		for (int ct = 0; ct <= n; ct++) {//A for loop to go through each fibonnaci Number till count is over
			if (ct <= 1) {
				cout << "Fibonnaci " << ct << " Result " << ct << endl;//If the number is less than or equal to 1 it will output 1
			}
			else {
				fib = fib1 + fib2;
				fib2 = fib1;
				fib1 = fib;
				cout << "Fibonnaci " << ct << " Result " << fib << endl;//Otherwise it performs fibonnaci and then outputs the result.
			}
		}
		return 0;
	}
	int FactorialRecursion(int n) {
		if (n > 1) return n * FactorialRecursion(n - 1);//Factorial Recrusion Code.
		else return 1;//Termination condition
	}
	int FactorialIterative(int n) {
		int res=1;
		for (int ct = 1; ct <= n; ct++) {//Goes through each factorial number till count is over
			res = res * ct;
			cout << "Factorial Of " << ct << " Is " << res << endl;//Returns the factorial using a totalling variable, res to output each.
		}
		return 0;
	}
	int UsingRecursion() {
		//The below utalised the functions above to time how long it takes for each code to run 
		//as well as looping the code for the recursion functions.
		clock_t t1 = clock();
		clock_t t2 = clock();
		int R1, R2, R3, R4;
		t1 = clock();
		cout << "Fibonacci Recurssion: " << endl;
		for (int ct = 0; ct <= 30; ct++)
			cout << "Fibonacci: " << ct << " is " << FibonacciRecursion(ct) << endl;

		t2 = clock();
		R1 = t2 - t1;
		cout << "Time taken is " << t2 - t1 << "\n";

		t1 = clock();

		cout << "\n";
		cout << "Factorial Recussion: " << endl;
		for (int ct = 0; ct <= 30; ct++)
			cout << "Factorial " << ct << " is " << FactorialRecursion(ct) << "\n";

		t2 = clock();
		R2 = t2 - t1;
		cout << "Time taken is " << t2 - t1 << "\n";

		t1 = clock();

		cout << "\n";
		cout << "Fibonacci Iterative" << endl;
		FibonacciIterative(30);

		t2 = clock();
		R3 = t2 - t1;
		cout << "Time taken is " << t2 - t1 << "\n";

		t1 = clock();

		cout << "\n";
		cout << "Factorial Iterative" << endl;
		FactorialIterative(30);

		t2 = clock();
		R4 = t2 - t1;
		cout << "Time taken is " << t2 - t1 << "\n";

		cout << "\n" << endl;
		cout << "Run Times: " << endl;
		cout << "Fibonacci Recurssion: " << R1 << endl;
		cout << "Factorial Recurssion: " << R2 << endl;
		cout << "Fibonacci Iterative: " << R3 << endl;
		cout << "Factorial Iterative: " << R4 << endl;
		return 0;
	}
}

int _tmain(int argc, _TCHAR* argv[]) {
	MatrixThings::UsingMatrix();
	RecursionWork::UsingRecursion();
	return 0;
}

// returns the factorial of integer n which is >= 0
// factorial of 0 or 1 is 1, otherwise it is n times factorial of n-1

namespace Week5Session {
	//NQueens
	//Hexagons
	//KnightsTour
	/*
	* Function tryRow(r)
	* //Try to put queen in row R
	*	for(int col=1;col<=n;col++)
	*		if(canputQueen(r,c)){
	*			putQueenThere(r,c)
	*			if(r<n tryRow(r+1)
	*				else solved = true
	*			if(Solved == false)
	*				removeQueenThere(r,c)
	*
	* Implementation
	* Board[r] = C
	* Specify that queen in column c in row R
	* canPutQueen(Row,Col)
	* How to test no Queen in this column
	*
	*/
}

/*Week 4 Tests
	//RecursionWork::FactorialIterative(10);
	RecursionWork::FibonacciInterative(10);
	//RecursionWork::FibonacciRecussion(10);
	//cout << "Factorial Example\nWritten by Richard Mitchell\n";
	// program calculates and outputs the factorial of 0 to 10
	//for (int ct = 0; ct <= 10; ct++)
	//	cout << "Factorial " << ct << " is " << RecursionWork::FactorialRecurssion(ct) << "\n";
	clock_t t1 = clock();
	//RecursionWork::OutputtingFibonacciRecussion();
	//RecursionWork::OutputtingFactorialRecurssion();
	clock_t t2 = clock();
	cout << "Time taken is " << t2 - t1 << "\n";

	int n = 9;

	cout << RecursionWork::fib(n);*/
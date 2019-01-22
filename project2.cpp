#include"stdafx.h"
#include <iostream>
#include <string>
#include <fstream>
#include<cstdlib>

#define MaxRows 500
#define MaxCols 500

using namespace std;

struct Image {
private:
	unsigned short SourceImage[MaxRows][MaxCols];
	int Rows, Cols;
	bool Available;

	void init() {
		Rows = 0;
		Cols = 0;
		Available = false;
	}
public:
	Image() {
		init();
	}

	Image(string Filename) {
		loadImage(Filename);
	}
	int SetPixel(unsigned short Value, int R, int C) {
		if (R > -1 && C > -1 && R < MaxRows && C < MaxCols) {
			SourceImage[R][C] = Value;
			return Value;
		}
		return -1;
	}

	bool isAvailable() {
		return Available;
	}
	string code, comment;
	int max;
	int loadImage(string FileName) {
		init();
		// Code to load an image from a storage into array named Image in main memory
		ifstream input(FileName);
		if (!input)
			return -1;

		getline(input, code);
		if (code != "P2")
			return -1;
		getline(input, comment);
		input >> Cols >> Rows >> max;
		for (int i = 0; i < Rows; i++)
			for (int j = 0; j < Cols; j++)
				input >> SourceImage[i][j];

		Available = true;
		return 0;
	}

	int saveImage(string FileName) {
		if (!Available)
			return 1;
		// Code to store the SourceImage or ProcessedImage into a file on some backing storage
		//  WhichImage 0 means copy source image and otherwise ProcessedImage
		ofstream output(FileName);

		output << code << endl;
		output << "# Created by Bajwa and Pal" << endl;
		output << Cols << " " << Rows << endl << max << endl;
		for (int i = 0; i < Rows; i++){
			for (int j = 0; j < Cols; j++)
				output << SourceImage[i][j] << " ";
			output << endl;
		}
		return 0;
	}

	void RotateImageAboutPoint(Image &Result, double ByAngle = 90, int Cx = 0, int Cy = 0) {
		// Code to compute a rotated version of original image
		// The image will be rotated about the center (Cx, Cy)
		// and the resultant image will be stored in Result
		int r, c;
		for (int i = 0; i < Rows; i++)
			for (int j = 0; j < Cols; j++)
				Result.SourceImage[i][j] = 0;

		double angle = (ByAngle*(22 / 7)) / 180;
		for (int i = 0; i < Rows; i++){
			for (int j = 0; j < Cols; j++){
				r = (cos(angle)*(i - Cx) - sin(angle)*(j - Cy)) + Cx;
				c = (sin(angle)*(i - Cx) + cos(angle)*(j - Cy)) + Cy;
				if (r >= 0 && c >= 0)
					Result.SourceImage[r][c] = SourceImage[i][j];
			}
		}
		Result.Available = true;
		Result.code = code;
		Result.comment = comment;
		Result.Cols = Cols;
		Result.Rows = Rows;
		Result.max = max;
	}

	void TranslateImage(Image &Result, int Tx, int Ty) {
		// Code to compute a translated version of original image
		// The image will be translated by Tx and Ty respectively along X and Y axes
		for (int i = 0; i < Rows; i++)
			for (int j = 0; j < Cols; j++)
				Result.SourceImage[i][j] = 0;

		for (int i = 0; i < Rows; i++)
			for (int j = 0; j < Cols; j++) 
				if ((i + Ty) >= 0 && (j + Tx) >= 0)
					Result.SourceImage[i + Ty][j + Tx] = SourceImage[i][j];
			
		Result.Available = true;
		Result.code = code;
		Result.comment = comment;
		Result.Cols = Cols + Tx;
		Result.Rows = Rows + Ty;
		Result.max = max;
	}

	void ScaleImage(Image &Result, double SX, double SY,int opt) {
		// Code to compute a scaled version of the original image
		// The image will be scaled by ratios Sx and Sy respectively along X and Y axes
		int r, c;
		if (opt == 2) {
			for (int i = 0, r = 0; i < Rows; i = i + SY, r++)
				for (int j = 0, c = 0; j < Cols; j = j + SX, c++)
					Result.SourceImage[r][c] = SourceImage[i][j];

			Result.Available = true;
			Result.code = code;
			Result.comment = comment;
			Result.Cols = Cols / SX;
			Result.Rows = Rows / SY;
			Result.max = max;
		}
		if (opt == 1)
		{
			for (double i = 0, r = 0; i < Rows; i = i + (1.0 / SY), r++)
				for (double j = 0, c = 0; j < Cols; j = j + (1.0 / SX), c++)
					Result.SourceImage[(int)r][(int)c] = SourceImage[(int)i][(int)j];

			Result.Available = true;
			Result.code = code;
			Result.comment = comment;
			Result.Cols = Cols * SX;
			Result.Rows = Rows * SY;
			Result.max = max;
		}
	}

	void LinearTransform(Image &Result, double TransformationMatrix[4]) {
		// Code to compute a transformed version of the original image
		// The transformation matrix will be specified in LinearTransformation matrix
		int r, c;
		for (int i = 0; i < Rows; i++)
			for (int j = 0; j < Cols; j++)
				Result.SourceImage[i][j] = 0;

		for (int i = 0; i < Rows; i++) {
			for (int j = 0; j < Cols; j++) {
				r = i*TransformationMatrix[0]+j*TransformationMatrix[1];
				c = i*TransformationMatrix[2]+j*TransformationMatrix[3];
				if (r >= 0 && c >= 0)
					Result.SourceImage[r][c] = SourceImage[i][j];
			}
		}
		Result.Available = true;
		Result.code = code;
		Result.comment = comment;
		Result.Cols = Cols;
		Result.Rows = Rows;
		Result.max = max;
	}

	void FlipImage(Image &Result, int Axis) {
		// Code to compute a flipped version of the original image
		// The Axis parameter will specify the flip axis.
		// 1 means flip along X-axis and otherwise flip along Y-axis
		if (Axis == 1)
			for (int i = 0; i < Rows; i++)
				for (int j = 0; j < Cols; j++)
					Result.SourceImage[i][Cols - j] = SourceImage[i][j];

		if (Axis == 2)
			for (int i = 0; i < Rows; i++)
				for (int j = 0; j < Cols; j++)
					Result.SourceImage[Rows - i][j] = SourceImage[i][j];

		Result.Available = true;
		Result.code = code;
		Result.comment = comment;
		Result.Cols = Cols;
		Result.Rows = Rows;
		Result.max = max;
	}

	void QuantizeImage(Image &Result, int RangeArray[]) {

		for (int i = 0; i < Rows; i++)
			for (int j = 0; j < Cols; j++)
				for (int k = 0; k < 8; k++){
					if ((RangeArray[k] - SourceImage[i][j]) >= 16 && (RangeArray[k] - SourceImage[i][j]) <= 32)
						Result.SourceImage[i][j] = RangeArray[k];
					if ((RangeArray[k] - SourceImage[i][j]) >= 0 && (RangeArray[k] - SourceImage[i][j]) <= 15)
						Result.SourceImage[i][j] = RangeArray[k - 1];
				}

		Result.Available = true;
		Result.code = code;
		Result.comment = comment;
		Result.Cols = Cols;
		Result.Rows = Rows;
		Result.max = max;

	}

	void MeanFilterImage(Image &Result) {
		for (int i = 0; i < Rows; i++){
			for (int j = 0; j < Cols; j++){
				if ((i + 1) >= Rows){
					SourceImage[i + 1][j] = 0;
					SourceImage[i + 1][j + 1] = 0;
					SourceImage[i + 1][j - 1] = 0;
				}
				if ((i - 1) < 0){
					SourceImage[i - 1][j] = 0;
					SourceImage[i - 1][j + 1] = 0;
					SourceImage[i - 1][j - 1] = 0;
				}
				if ((j + 1) >= Cols){
					SourceImage[i][j + 1] = 0;
					SourceImage[i + 1][j + 1] = 0;
					SourceImage[i - 1][j + 1] = 0;
				}
				if ((j - 1) < 0){
					SourceImage[i][j - 1] = 0;
					SourceImage[i - 1][j - 1] = 0;
					SourceImage[i + 1][j - 1] = 0;
				}
				Result.SourceImage[i][j] = (SourceImage[i - 1][j - 1] + SourceImage[i - 1][j] +
					SourceImage[i - 1][j + 1] + SourceImage[i][j - 1] + SourceImage[i][j] +
					SourceImage[i][j + 1] + SourceImage[i + 1][j - 1] + SourceImage[i + 1][j] +
					SourceImage[i + 1][j + 1]) / 9;
			}
		}
		Result.Available = true;
		Result.code = code;
		Result.comment = comment;
		Result.Cols = Cols;
		Result.Rows = Rows;
		Result.max = max;
	}

	void MedianFilterImage(Image &Result) {
		for (int i = 0; i < Rows; i++)
			for (int j = 0; j < Cols; j++)
			{
				if ((i + 1) >= Rows)
				{
					SourceImage[i + 1][j] = 0;
					SourceImage[i + 1][j + 1] = 0;
					SourceImage[i + 1][j - 1] = 0;
				}
				if ((i - 1) < 0)
				{
					SourceImage[i - 1][j] = 0;
					SourceImage[i - 1][j + 1] = 0;
					SourceImage[i - 1][j - 1] = 0;
				}
				if ((j + 1) >= Cols)
				{
					SourceImage[i][j + 1] = 0;
					SourceImage[i + 1][j + 1] = 0;
					SourceImage[i - 1][j + 1] = 0;
				}
				if ((j - 1) < 0)
				{
					SourceImage[i][j - 1] = 0;
					SourceImage[i - 1][j - 1] = 0;
					SourceImage[i + 1][j - 1] = 0;
				}
				int Array[9];
				Array[0] = SourceImage[i - 1][j];
				Array[1] = SourceImage[i - 1][j - 1];
				Array[2] = SourceImage[i - 1][j + 1];
				Array[3] = SourceImage[i][j];
				Array[4] = SourceImage[i][j - 1];
				Array[5] = SourceImage[i][j + 1];
				Array[6] = SourceImage[i + 1][j];
				Array[7] = SourceImage[i + 1][j - 1];
				Array[8] = SourceImage[i + 1][j + 1];

				for (int k = 0; k < 9; k++)
				{
					int min = k;
					for (int m = k + 1; m < 9; m++)
						if (Array[min] > Array[m])
							min = m;
					if (Array[min] != Array[k])
						swap(Array[min], Array[k]);
				}
				Result.SourceImage[i][j] = Array[4];
			}

		Result.Available = true;
		Result.code = code;
		Result.comment = comment;
		Result.Cols = Cols;
		Result.Rows = Rows;
		Result.max = max;


	}

	void ComputeDerivativeImage(Image &Result) {
		for (int i = 0; i < Rows; i++) {
			for (int j = 0; j < Cols; j++) {
				if ((i + 1) >= Rows){
					SourceImage[i + 1][j] = 0;
					SourceImage[i + 1][j + 1] = 0;
					SourceImage[i + 1][j - 1] = 0;
				}
				if ((i - 1) < 0){
					SourceImage[i - 1][j] = 0;
					SourceImage[i - 1][j + 1] = 0;
					SourceImage[i - 1][j - 1] = 0;
				}
				if ((j + 1) >= Cols){
					SourceImage[i][j + 1] = 0;
					SourceImage[i + 1][j + 1] = 0;
					SourceImage[i - 1][j + 1] = 0;
				}
				if ((j - 1) < 0){
					SourceImage[i][j - 1] = 0;
					SourceImage[i - 1][j - 1] = 0;
					SourceImage[i + 1][j - 1] = 0;
				}
				int tmp = SourceImage[i - 1][j - 1] * (1) +
					SourceImage[i - 1][j] * (1) + SourceImage[i - 1][j + 1] * (1)
					+ SourceImage[i][j - 1] * (1) + SourceImage[i][j + 1] * (-1)
					+ SourceImage[i + 1][j - 1] * (-1) + SourceImage[i + 1][j] * (-1)+
					SourceImage[i + 1][j + 1] * (-1);
				if (tmp < 0)
					tmp = -1 * tmp;
				Result.SourceImage[i][j] =tmp ;
			}
		}
		Result.Available = true;
		Result.code = code;
		Result.comment = comment;
		Result.Cols = Cols;
		Result.Rows = Rows;
		Result.max = max;
	}

	void ApplyFilter(Image &Result) {
		for (int i = 0; i < Rows; i++) {
			for (int j = 0; j < Cols; j++) {
				if ((i + 1) >= Rows) {
					SourceImage[i + 1][j] = 0;
					SourceImage[i + 1][j + 1] = 0;
					SourceImage[i + 1][j - 1] = 0;
				}
				if ((i - 1) < 0) {
					SourceImage[i - 1][j] = 0;
					SourceImage[i - 1][j + 1] = 0;
					SourceImage[i - 1][j - 1] = 0;
				}
				if ((j + 1) >= Cols) {
					SourceImage[i][j + 1] = 0;
					SourceImage[i + 1][j + 1] = 0;
					SourceImage[i - 1][j + 1] = 0;
				}
				if ((j - 1) < 0) {
					SourceImage[i][j - 1] = 0;
					SourceImage[i - 1][j - 1] = 0;
					SourceImage[i + 1][j - 1] = 0;
				}
				double sum1, sum2;
				int sum;
				sum1 = (SourceImage[i - 1][j - 1] + SourceImage[i - 1][j]
					+ SourceImage[i - 1][j + 1] + SourceImage[i][j - 1]
					+ SourceImage[i][j + 1] + SourceImage[i + 1][j - 1]
					+ SourceImage[i + 1][j] + SourceImage[i + 1][j + 1])*(0.111);
				sum2 = SourceImage[i][j] * (1.888);
				sum = sum1 - sum2;
				if (sum < 0)
					sum = sum*-1;
				Result.SourceImage[i][j] = sum ;
			}
		}
		Result.Available = true;
		Result.code = code;
		Result.comment = comment;
		Result.Cols = Cols;
		Result.Rows = Rows;
		Result.max = max;
	}
};

struct Menu {
private:
	string Options[30];
	int MaxOptions;
public:

	int GetMaxOptions() {
		return MaxOptions;
	}

	Menu(string FileName) {
		MaxOptions = 0;

		ifstream OptionsFile(FileName.c_str());

		if (!OptionsFile) {
			return;
		}
		while (!OptionsFile.eof()) {
			getline(OptionsFile, Options[MaxOptions++]);
		}
		return;
	}

	int ShowMenuAndGetChoice() {
		int Choice = 1;

		cout << "\nWel-Come To FUN Mini-Project About Image Processing (Arrays)\n" <<
			"Created by Bajwa and Pal.\n";
		do
		{
			if (Choice < 1 || Choice > MaxOptions) {
				cout << endl << " Please Make a Valid Choice\n\n Press enter to continue";
			}
			for (int i = 0; i< MaxOptions; i++)
				cout << endl << i + 1 << ":\t" << Options[i];

			cout << endl << endl << MaxOptions + 1 << ":\t Exit";

			cout << endl << "\nMake a choice from the Menu by Specifying Choice No ";
			cin >> Choice;
		} while (Choice < 1 || Choice > MaxOptions + 1);

		return Choice;
	}
};

void LOAD(Image &I) {
	string FileName = "";
	cout << "Specify Image File Name: ";
	do {
		getline(cin, FileName);
	} while (FileName.length() < 1);

	if (I.loadImage(FileName) != 0) {
		cout << endl << "File Error: Image Not Loaded" << endl << endl;
	}
	else {
		cout << endl << " Image has Been Loaded" << endl << endl;
	}
}

void SAVE(Image &I) {
	if (!I.isAvailable()) {
		cout << endl << "No Image Available" << endl << endl;
		return;
	}
	string FileName = "";
	cout << "Specify Image File Name: ";
	do {
		getline(cin, FileName);
	} while (FileName.length() < 1);

	int Result = I.saveImage(FileName);

	if (2 == Result) {
		cout << endl << "File Error: Image Not Saved" << endl << endl;
	}
	else {
		cout << endl << " Image Saved Successfully" << endl << endl;
	}
	return;
}

void Rotate(Image &I, Image &Result) {
	// Ask user to specify starting and ending rotation angle and a step size
	// compute and save in different files all rotated versions of the source image I.
	int Cx, Cy;
	cout << "enter center co-ordinates (Cx,Cy):";
	cin >> Cx >> Cy;
	int angle;
	cout << "enter angle in degrees:";
	cin >> angle;
	I.RotateImageAboutPoint(Result, angle, Cx, Cy);
}

void Translate(Image &I, Image &Result) {
	// Ask user to specify Translation in X and Y and a step size
	// compute and save translated images in different files that are translated
	// versions of the source image I.
	int Tx, Ty;
	cout << "enter translation(X,Y):";
	cin >> Tx >> Ty;
	I.TranslateImage(Result, Tx, Ty);
}

void Flip(Image &I, Image &Result) {
	// Main Driver for computing Flipped version of image I
	// It must save the flipped version into image Result
	int axis;
	cout << "enter axis:\n" << "1. Along X axis\n" << "2. Along Y axis\n";
	cin >> axis;
	I.FlipImage(Result, axis);
}

void Scale(Image &I, Image &Result) {
	// Main Driver for computing Flipped version of image I
	// It must save the flipped version into image Result
	double Sx, Sy;
	int opt;
	cout << "enter choice:\n" << "1.zoom in\n 2.zoom put\n";
	cin >> opt;
	cout << "enter scale ratio (Sx,Sy):";
	cin >> Sx >> Sy;
	I.ScaleImage(Result, Sx, Sy,opt);
}

void Transform(Image &I, Image &Result) {
	double TransformationMatrix[4];
	cout << "enter 4 values for transformation:";
	for (int i = 0; i < 4; i++)
		cin >> TransformationMatrix[i];
	I.LinearTransform(Result,TransformationMatrix);
}

void Quantize(Image &I, Image &Result) {
	// Main Driver for Quantizing an Image
	// It must compute Quantized image by quantizing image I and save result in image Result
	int RangeArray[9];
	int k = 0;
	for (int i = 0; i <= 255; i++){
		if (k <= 255)
			RangeArray[i] = k;
		k = k + 32;
	}
	I.QuantizeImage(Result, RangeArray);
}

void MeanFilter(Image &I, Image &Result) {
	// Main driver for mean filter
	// It must apply a 3 x 3 mean filter on image I and save result in image Result
	I.MeanFilterImage(Result);
}

void MedianFilter(Image &I, Image &Result) {
	// Main Driver for Median  Filter
	// It must apply 3 x 3 median filter on image I and Save Result in Image Result
	I.MedianFilterImage(Result);
}

void ComputeDerivative(Image &I, Image &Result) {
	// Main Driver for Computing Derivative
	// It must Compute Derivative of image I and Save Result in Image Result
	I.ComputeDerivativeImage(Result);
}

void Filter(Image &I, Image &Result) {
	// Main Driver for Filter
	// It must apply a filter on image I and Save Result in Image Result
	I.ApplyFilter(Result);
}

int main() {
	Image A, B;
	int Choice = -2;
	Menu MainMenu("MainMenu.txt");

	do {
		Choice = MainMenu.ShowMenuAndGetChoice();

		if (1 == Choice) {
			LOAD(A);
		}
		else if (2 == Choice) {
			SAVE(A);
		}
		else if (3 == Choice) {
			SAVE(B);
		}
		else if (4 == Choice) {
			Rotate(A, B);
		}
		else if (5 == Choice) {
			Translate(A, B);
		}
		else if (6 == Choice) {
			Scale(A, B);
		}
		else if (7 == Choice) {
			Transform(A, B);
		}
		else if (8 == Choice) {
			Flip(A, B);
		}
		else if (9 == Choice) {
			Quantize(A, B);
		}
		else if (10 == Choice) {
			MeanFilter(A, B);
		}
		else if (11 == Choice) {
			MedianFilter(A, B);
		}
		else if (12 == Choice) {
			ComputeDerivative(A, B);
		}
		else if (13 == Choice) {
			Filter(A, B);
		}
		else{
			cout << endl << "\n\tLeaving already. \n\tHave FUN. \n\tSee You Soon\t";
			exit(0);
		}
	} while (Choice != (MainMenu.GetMaxOptions() + 1));
	return 0;
}


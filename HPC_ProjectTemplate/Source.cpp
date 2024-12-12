#include <iostream>
#include <math.h>
#include <stdlib.h>
#include<string.h>
#include<msclr\marshal_cppstd.h>
#include <ctime>// include this header 
#include <vector>
#include <cmath>
#include <limits>
#include <cstdlib>
#pragma once

#using <mscorlib.dll>
#using <System.dll>
#using <System.Drawing.dll>
#using <System.Windows.Forms.dll>
#include <mpi.h>
using namespace std;
using namespace msclr::interop;
using namespace System;
using namespace System::Windows::Forms;

int* inputImage(int* w, int* h, System::String^ imagePath) //put the size of image in w & h
{
	int* input;
	int OriginalImageWidth, OriginalImageHeight;
	//*********************************************************Read Image and save it to local arrayss*************************	
	//Read Image and save it to local arrayss

	System::Drawing::Bitmap BM(imagePath);

	OriginalImageWidth = BM.Width;
	OriginalImageHeight = BM.Height;
	*w = BM.Width;
	*h = BM.Height;
	int *Red = new int[BM.Height * BM.Width];
	int *Green = new int[BM.Height * BM.Width];
	int *Blue = new int[BM.Height * BM.Width];
	input = new int[BM.Height*BM.Width];
	for (int i = 0; i < BM.Height; i++)
	{
		for (int j = 0; j < BM.Width; j++)
		{
			System::Drawing::Color c = BM.GetPixel(j, i);

			Red[i * BM.Width + j] = c.R;
			Blue[i * BM.Width + j] = c.B;
			Green[i * BM.Width + j] = c.G;

			input[i*BM.Width + j] = ((c.R + c.B + c.G) / 3); //gray scale value equals the average of RGB values

		}

	}
	return input;
}


void createImage(int* image, int width, int height, int index)
{
	System::Drawing::Bitmap MyNewImage(width, height);


	for (int i = 0; i < MyNewImage.Height; i++)
	{
		for (int j = 0; j < MyNewImage.Width; j++)
		{
			//i * OriginalImageWidth + j
			if (image[i*width + j] < 0)
			{
				image[i*width + j] = 0;
			}
			if (image[i*width + j] > 255)
			{
				image[i*width + j] = 255;
			}
			System::Drawing::Color c = System::Drawing::Color::FromArgb(image[i*MyNewImage.Width + j], image[i*MyNewImage.Width + j], image[i*MyNewImage.Width + j]);
			MyNewImage.SetPixel(j, i, c);
		}
	}
	MyNewImage.Save("..//Data//Output//outputRes" + index + ".png");
	std::cout << "result Image Saved " << index << endl;
}

void displayImage(System::String^ imagePath) {
	// Create a new Windows Form
	Form^ form = gcnew Form();
	form->Text = "Input Image";
	form->Size = System::Drawing::Size(800, 600);

	// Create a PictureBox
	PictureBox^ pictureBox = gcnew PictureBox();
	pictureBox->Dock = DockStyle::Fill;
	pictureBox->SizeMode = PictureBoxSizeMode::Zoom;

	// Load the image into the PictureBox
	pictureBox->Image = System::Drawing::Image::FromFile(imagePath);

	// Add PictureBox to the form
	form->Controls->Add(pictureBox);

	// Show the form
	Application::Run(form);
}

int calcNumberOfchunkSize(int imglength,int size) {
	//int chunks = imglength / size;
	//if (chunks % size != 0) {
	//	chunks++;
	//}
	return  imglength / size;
}

//(x,y) (x2,y2)
int euclideanDistance(int point1,int point2) {
	return sqrt((point1 - point2) * (point1 - point2));
}

int main()
{
	int ImageWidth = 4, ImageHeight = 4;
	int start_s, stop_s, TotalTime = 0;
	System::String^ imagePath;
	std::string img;
	img = "..//Data//Input//test.png";
	imagePath = marshal_as<System::String^>(img);
	int* imageData = inputImage(&ImageWidth, &ImageHeight, imagePath); //array of intinsities
	start_s = clock();
	/***********************************************************/

	MPI_Init(NULL, NULL);
	int size; //number of processors
	int rank; //id of each processor
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int k = 3;
	int chunkSize = calcNumberOfchunkSize(ImageHeight * ImageWidth, size);//number assigned to processors
	int* recvchunk = new int[chunkSize];
	int* centroids= new int[3]();
	int* localdata= new int [chunkSize];
	int* newCentroids= new int[3]();
	bool converge = false;
	int itr = 100;
	int *count = new int[3]();
	int* sum = new int[3]();
	int* localSum = new int[3]();
	int* localCount = new int[3]();
	int* globalImage=new int[size * chunkSize];
	MPI_Scatter(imageData, chunkSize, MPI_INT, recvchunk, chunkSize, MPI_INT, 0, MPI_COMM_WORLD);
	
	if (rank == 0) {
		for (size_t i = 0; i < 3; i++)
		{
			centroids[i] = (i*255)/(2);
		}
	}
	MPI_Bcast(centroids, 3, MPI_INT, 0, MPI_COMM_WORLD);

	while (!converge&&itr>0) { //while(true)
		fill(localSum, localSum + 3, 0);
		fill(localCount, localCount + 3, 0);

		for (int i = 0; i < chunkSize; ++i) { 
			int minDistance = 1000000;
			int nearest = 0;
			for (int j = 0; j < 3; ++j) {
				int distance = euclideanDistance(recvchunk[i], centroids[j]);
				if (distance < minDistance) {
					minDistance = distance;
					nearest = j;
				}
				localdata[i] = centroids[nearest];
				localCount[nearest]++;
				localSum[nearest] += recvchunk[i];
			}
		}
		MPI_Reduce(localSum, sum, 3, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(localCount, count, 3, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

		if (rank == 0) {
		
			for (size_t i = 0; i < 3; i++)
			{
				newCentroids[i] = sum[i] / count[i];
				if (abs(newCentroids[i] - centroids[i]) > 0.5) {
					converge = true;
				}
			}
			for (size_t i = 0; i < 3; i++)
			{
				centroids[i] = newCentroids[i];
			}
		}
		MPI_Bcast(centroids, k, MPI_INT, 0, MPI_COMM_WORLD);
		itr--;
	}
	MPI_Gather(localdata, chunkSize, MPI_INT, globalImage, chunkSize, MPI_INT, 0, MPI_COMM_WORLD);
	if(rank ==0)
	{
		for (int i = ((ImageHeight * ImageWidth) / size) * size; i < ImageHeight * ImageWidth; ++i) { //ely gowa recvchunk
			int minDistance = 1000000;
			int nearest = 0;
			for (int j = 0; j < 3; ++j) {
				int distance = euclideanDistance(imageData[i], centroids[j]);
				if (distance < minDistance) {
					minDistance = distance;
					nearest = j;
				}
				localdata[i] = centroids[nearest];
			}
			imageData[i] = centroids[nearest];

		}
		for (size_t i = 0; i < 3; i++)
		{
			cout << centroids[i] << " ";
		}
		cout << ((ImageHeight * ImageWidth) / size) * size <<endl;
		cout << (ImageHeight * ImageWidth);
		for (int i = ((ImageHeight * ImageWidth) / size) * size; i < ImageHeight * ImageWidth; ++i) {
			cout << imageData[i] << " ";
		}
		for (size_t i = 0; i < chunkSize; i++)
		{
			cout << localdata[i]<<" ";
		}
	}

	if (rank == 0) {
		createImage(globalImage, ImageWidth, ImageHeight, 1);
	}
	delete[] imageData;
	delete[] recvchunk;
	delete[] localSum;
	delete[] localCount;
	delete[] sum;
	delete[] count;

	MPI_Finalize();

	/***********************************************************/
	stop_s = clock();
	TotalTime += (stop_s - start_s) / double(CLOCKS_PER_SEC) * 1000;
	std::cout << "time: " << TotalTime << endl;

	return 0;

}
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

int calcNumberOfChunks(int imglength,int size) {
	int chunks = imglength / size;
	if (chunks % 2 != 0) {
		chunks++;
	}
	return chunks;
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
	//displayImage(imagePath);
	int k = 3;
	int chunks = calcNumberOfChunks(ImageHeight * ImageWidth, size);//number assigned to processors
	int* recvchunk = new int[chunks];
	int w=0, g=0, b=0;
	int centroids[3] = { 255, 125, 0};// Initial centroids
	int* labels = new int[ImageHeight * ImageWidth];
	MPI_Scatter(imageData, chunks, MPI_INT, recvchunk, chunks, MPI_INT, 0, MPI_COMM_WORLD);

	while(true)
	{
		vector<int> whiteArr, grayArr, blackArr;
		for (int i = 0; i < chunks; ++i) {
			int minDistance = 300;
			for (int j = 0; j < k; ++j) {
				int distance = euclideanDistance(recvchunk[i], centroids[j]);
				if (distance < minDistance) {
					minDistance = distance;
					labels[i + rank * chunks] = j;
					if (labels[i + rank * chunks] == 0) whiteArr.push_back(recvchunk[i]);
					else if (labels[i + rank * chunks] == 1) grayArr.push_back(recvchunk[i]);
					else blackArr.push_back(recvchunk[i]);

				}
			}
		}
		if (rank == 0) {
			int t = 2;
			while (t >= 0)
			{
				int sum = 0;
				if (t == 0)
				{
					if(whiteArr.size()>0)
					{
						for (size_t i = 0; i < whiteArr.size(); i++)
						{
							sum += whiteArr[i];
						}
						w = centroids[0];
						centroids[0] = sum / (whiteArr.size());
					}
				}
				else if (t == 1)
				{
					if(grayArr.size()>0)
					{
						for (size_t i = 0; i < grayArr.size(); i++)
						{
							sum += grayArr[i];
						}
						g = centroids[1];
						centroids[1] = sum / (grayArr.size());
					}
				}
				else
				{
					if(blackArr.size()>0)
					{
						for (size_t i = 0; i < blackArr.size(); i++)
						{
							sum += blackArr[i];
						}
						b = centroids[2];
						centroids[2] = sum / (blackArr.size());
					}
				}
				t--;
			}
			if (centroids[0] == w && centroids[1] == g && centroids[2] == b) {	
				break;
			}
			MPI_Bcast(centroids, k, MPI_INT, 0, MPI_COMM_WORLD);
		}
	}
	MPI_Gather(labels, chunks, MPI_INT, imageData, chunks, MPI_INT, 0, MPI_COMM_WORLD);
	if(rank==0)
	{
		for (int i = 0; i < ImageHeight * ImageWidth; i++)
		{
			if (labels[i] == 0)labels[i] = centroids[0];
			else if(labels[i]==1)labels[i] = centroids[1];
			else labels[i] = centroids[2];
		}
	}
	MPI_Finalize();
	/***********************************************************/
	stop_s = clock();
	TotalTime += (stop_s - start_s) / double(CLOCKS_PER_SEC) * 1000;
	createImage(labels, ImageWidth, ImageHeight, 1);
	std::cout << "time: " << TotalTime << endl;
	free(imageData);
	return 0;

}




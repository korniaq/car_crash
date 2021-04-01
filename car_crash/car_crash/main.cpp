#include"mesh.h"
#include"solve.h"
#define _CRT_SECURE_NO_WARNINGS

int main()
{
	graphics(1200, 600);
	scale(-400, -100, 600, 400);
	Mesh<double> mesh(10);
	mesh.ReadNodes("car.txt");
	mesh.CreateMesh();

	const int iter = 20;
	wait();
	//Solve(mesh, iter);
	Animate(mesh, iter);
	wait();

	return 0;
}
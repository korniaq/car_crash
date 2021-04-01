#pragma once

#include"winbgi2.h"
#include<iostream>
#include<vector>
#include<algorithm>
#include<cmath>
#include<fstream>
#include<array>

///////////////////////////////////////////
//class Node
template<class T>
class Node
{
private:
	T x;
	T y;
	int id;
	bool fix;
	bool ground;
	template<class T> friend class Element;
	template<class T> friend class Mesh;
public:
	///////////////////////////////////////////
	//creating stuff
	Node(const T& xx = 0, const T& yy = 0, const int& idid = 0, const bool& fixx = 0, const bool& grnd = 0) :
		x(xx), y(yy), id(idid), fix(fixx), ground(grnd) {}
	~Node() = default;

	///////////////////////////////////////////
	//calculating stuff
	const T& CalculateDistance(const T& xx, const T& yy)
	{
		return sqrt((x - xx) * (x - xx) + (y - yy) * (y - yy));
	}

	///////////////////////////////////////////
	//operators
	Node<T>& operator=(const Node<T>&) = default; 
	bool operator==(const Node<T>& theother_node)
	{
		if (this->x == theother_node.x && this->y == theother_node.y)
			return true;
		else return false;
	}
	template<class T> friend std::ostream& operator<<(std::ostream&, const Node<T>&);
};

template<class T>
std::ostream& operator<<(std::ostream& o, const Node<T>& node)
{
	return o << "id: " << node.id << "\t" << node.x << "\t" << node.y << "\t" << node.fix << std::endl;
}

///////////////////////////////////////////
//class Element
template<class T>
class Element
{
private:
	static constexpr double E = 2 * 10e6;
	static constexpr double A = 0.08;
	static constexpr double v = 0.003;
	static constexpr double density = 7.860;
	int id;
	std::array<Node<T>,3> nodes;
	T Area;
	T Mass;
	T Stress;
	template<class T> friend class Mesh;
public:
	std::array<std::array<T, 6>, 6> local_k;
	std::array<std::array<T, 6>, 6> local_m;

	///////////////////////////////////////////
	//creating stuff
	Element() = default;
	Element(const Node<T> & node0, const Node<T> node1, const Node<T> & node2, const int& idid = 0) : id(idid), Area(0), Mass(0), Stress(0)
	{
		nodes[0] = node0;
		nodes[1] = node1;
		nodes[2] = node2;
	}
	~Element() = default;

	///////////////////////////////////////////
	//calculating stuff
	const T& CalculateLength(const int& node1, const int& node2)
	{
		return sqrt((nodes[node1].x - nodes[node2].x) * (nodes[node1].x - nodes[node2].x) + (nodes[node1].y - nodes[node2].y) * (nodes[node1].y - nodes[node2].y));
	}
	Node<T>& CalculateMidNode(const int& node1, const int& node2, const int& no)
	{
		bool fix = false;
		if (nodes[node1].fix == true && nodes[node2].fix == true)
			fix = true;
		return *new Node<T>((nodes[node1].x + nodes[node2].x) / 2, (nodes[node1].y + nodes[node2].y) / 2, no, fix);
	}
	const T& TheLongestSide(int& node0, int& node1, int& node2)
	{
		T l1 = CalculateLength(0, 1);
		T l2 = CalculateLength(1, 2);
		T l3 = CalculateLength(2, 0);
		if (l1 > l2)
			if (l1 > l3)
			{ node0 = 0; node1 = 1; node2 = 2; return l1; }
			else
			{ node0 = 2; node1 = 0; node2 = 1; return l3; }
		else if (l2 > l3)
		{ node0 = 1; node1 = 2; node2 = 0; return l2; }
		else
		{ node0 = 2; node1 = 0; node2 = 1; return l3; }
	}
	void CalculateProperties()
	{
		T l1 = CalculateLength(0, 1);
		T l2 = CalculateLength(1, 2);
		T l3 = CalculateLength(2, 0);

		//area and mass
		T p = 0.5 * (l1 + l2 + l3);
		Area = sqrt(p * (p - l1) * (p - l2) * (p - l3));
		Mass = Area * density;

		//local stiffness matrix
		CalculateLocalStiffnessMatrix();

		//local mass matrix
		CalculateLocalMassMatrix();
	}
	void CalculateLocalStiffnessMatrix()
	{
		T l1 = CalculateLength(0, 1);
		T l2 = CalculateLength(1, 2);
		T l3 = CalculateLength(2, 0);

		T c1 = (nodes[1].x - nodes[0].x) / l1;
		T c2 = (nodes[2].x - nodes[1].x) / l2;
		T c3 = (nodes[0].x - nodes[2].x) / l3;
		T s1 = (nodes[1].y - nodes[0].y) / l1;
		T s2 = (nodes[2].y - nodes[1].y) / l2;
		T s3 = (nodes[0].y - nodes[2].y) / l3;

		local_k[0][0] = c1 * c1 / l1 + c3 * c3 / l3;
		local_k[0][1] = s1 * c1 / l1 + s3 * c3 / l3;
		local_k[0][2] = -c1 * c1 / l1;
		local_k[0][3] = -s1 * c1 / l1;
		local_k[0][4] = -c3 * c3 / l3;
		local_k[0][5] = -s3 * c3 / l3;

		local_k[1][0] = s1 * c1 / l1 + s3 * c3 / l3;
		local_k[1][1] = s1 * s1 / l1 + s3 * s3 / l3;
		local_k[1][2] = -s1 * c1 / l1;
		local_k[1][3] = -s1 * s1 / l1;
		local_k[1][4] = -s3 * c3 / l3;
		local_k[1][5] = -s3 * s3 / l3;

		local_k[2][0] = -c1 * c1 / l1;
		local_k[2][1] = -s1 * c1 / l1;
		local_k[2][2] = c1 * c1 / l1 + c2 * c2 / l2;
		local_k[2][3] = s1 * c1 / l1 + s2 * c2 / l2;
		local_k[2][4] = -c2 * c2 / l2;
		local_k[2][5] = -s2 * c2 / l2;

		local_k[3][0] = -s1 * c1 / l1;
		local_k[3][1] = -s1 * s1 / l1;
		local_k[3][2] = s1 * c1 / l1 + s2 * c2 / l2;
		local_k[3][3] = s1 * s1 / l1 + s2 * s2 / l2;
		local_k[3][4] = -s2 * c2 / l2;
		local_k[3][5] = -s2 * s2 / l2;

		local_k[4][0] = -c3 * c3 / l3;
		local_k[4][1] = -s3 * c3 / l3;
		local_k[4][2] = -c2 * c2 / l2;
		local_k[4][3] = -s2 * c2 / l2;
		local_k[4][4] = c2 * c2 / l2 + c3 * c3 / l3;
		local_k[4][5] = s2 * c2 / l2 + s3 * c3 / l3;

		local_k[5][0] = -s3 * c3 / l3;
		local_k[5][1] = -s3 * s3 / l3;
		local_k[5][2] = -s2 * c2 / l2;
		local_k[5][3] = -s2 * s2 / l2;
		local_k[5][4] = s2 * c2 / l2 + s3 * c3 / l3;
		local_k[5][5] = s2 * s2 / l2 + s3 * s3 / l3;

		for (int i = 0; i < 6; i++)
			for (int j = 0; j < 6; j++)
				local_k[i][j] *= E * A;

		//boundary conditions
		for (int i = 0; i < 3; i++)
			if (nodes[i].fix == true)
			{
				for (int j = 0; j < 6; j++)
				{
					local_k[2 * i][j] = 0;
					local_k[2 * i + 1][j] = 0;
					local_k[j][2 * i] = 0;
					local_k[j][2 * i + 1] = 0;
				}
				local_k[2 * i][2 * i] = 1;
				local_k[2 * i + 1][2 * i + 1] = 1;
			}
	}
	void CalculateLocalMassMatrix()
	{	
		T x0, y0, x1, y1, x2, y2;
		x0 = y0 = 0;
		x1 = nodes[1].x - nodes[0].x;
		y1 = nodes[1].y - nodes[0].y;
		x2 = nodes[2].x - nodes[0].x;
		y2 = nodes[2].y - nodes[0].y;
		T J = x0 * y1 + x1 * y2 + x2 * y0 - x2 * y1 - x1 * y0 - x0 * y2;

		T Nc[6] = { x1 * y2 - x2 * y1, x1 * y2 - x2 * y1, x2 * y0 - x0 * y2, x2 * y0 - x0 * y2, x0 * y1 - x1 * y0, x0 * y1 - x1 * y0 };
		T Nx[6] = { y1 - y2, y1 - y2, y2 - y0, y2 - y0, y0 - y1, y0 - y1 };
		T Ny[6] = { x2 - x1, x2 - x1, x0 - x2, x0 - x2, x1 - x0, x1 - x0 };

		//local_mij = ro * integral(Ni * Nj dxdy), where Ni = (Nci + Nxi*x + Nyi*y) / J
		//integrating with Gauss method, assumed x0,y0 = 0 for calculations to be less complicated
		for(int i = 0; i < 6; i ++)
			for (int j = 0; j < 6; j++)
			{
				local_m[i][j] = 0;
				local_m[i][j] += Nc[i] * Nc[j];
				local_m[i][j] += 1. / 3 * (Nc[i] * Nx[j] + Nx[i] * Nc[j]) * (x1 + x2);
				local_m[i][j] += 1. / 3 * (Nc[i] * Ny[j] + Ny[i] * Nc[j]) * (y1 + y2);
				local_m[i][j] += 1. / 6 * Nx[i] * Nx[j] * (x1 * x1 + x1 * x2 + x2 * x2);
				local_m[i][j] += 1. / 6 * Ny[i] * Ny[j] * (y1 * y1 + y1 * y2 + y2 * y2);
				local_m[i][j] += 1. / 6 * (Nx[i] * Ny[j] + Ny[i] * Nx[j]) * (x1 * y1 + 0.5 * x1 * y2 + 0.5 * y1 * x2 + x2 * y2);
				local_m[i][j] = local_m[i][j] * density / (2 * J);
			}

		//boundary conditions
		for (int i = 0; i < 3; i++)
			if (nodes[i].fix == true)
			{
				for (int j = 0; j < 6; j++)
				{
					local_m[2 * i][j] = 0;
					local_m[2 * i + 1][j] = 0;
					local_m[j][2 * i] = 0;
					local_m[j][2 * i + 1] = 0;
				}
				local_m[2 * i][2 * i] = 1;
				local_m[2 * i + 1][2 * i + 1] = 1;
			}
	}
	void CalculateStress(const T& dx0, const T& dy0, const T& dx1, const T& dy1, const T& dx2, const T& dy2)
	{
		T epsilonx, epsilony, gammaxy, sigmax, sigmay, sigmaz, tauxy;
		T x0, y0, x1, y1, x2, y2;
		x0 = y0 = 0;
		x1 = nodes[1].x - nodes[0].x;
		y1 = nodes[1].y - nodes[0].y;
		x2 = nodes[2].x - nodes[0].x;
		y2 = nodes[2].y - nodes[0].y;

		T ux0, uy0, ux1, uy1, ux2, uy2;
		ux0 = uy0 = 0;
		ux1 = dx1 - dx0;
		uy1 = dy1 - dy0;
		ux2 = dx2 - dx0;
		uy2 = dy2 - dy0;

		T dN0x, dN0y, dN1x, dN1y, dN2x, dN2y;
		dN0x = (y1 - y2) / (x0 * y1 + x1 * y2 + x2 * y0 - x2 * y1 - x1 * y0 - x0 * y2);
		dN0y = (x2 - x1) / (x0 * y1 + x1 * y2 + x2 * y0 - x2 * y1 - x1 * y0 - x0 * y2);
		dN1x = (y2 - y0) / (x0 * y1 + x1 * y2 + x2 * y0 - x2 * y1 - x1 * y0 - x0 * y2);
		dN1y = (x0 - x2) / (x0 * y1 + x1 * y2 + x2 * y0 - x2 * y1 - x1 * y0 - x0 * y2);
		dN2x = (y0 - y1) / (x0 * y1 + x1 * y2 + x2 * y0 - x2 * y1 - x1 * y0 - x0 * y2);
		dN2y = (x1 - x0) / (x0 * y1 + x1 * y2 + x2 * y0 - x2 * y1 - x1 * y0 - x0 * y2);

		epsilonx = dN0x * ux0 + dN1x * ux1 + dN2x * ux2;
		epsilony = dN0y * uy0 + dN1y * uy1 + dN2y * uy2;
		gammaxy = dN0x * uy0 + dN0y * ux0 + dN1x * uy1 + dN1y * ux1 + dN2x * uy2 + dN2y * ux2;

		sigmax = E / (1 + v) * (epsilonx + v / (1 - 2 * v) * (epsilonx + epsilony));
		sigmay = E / (1 + v) * (epsilony + v / (1 - 2 * v) * (epsilonx + epsilony));
		sigmaz = E * v / ((1 + v) * (1 - 2 * v)) * (epsilonx + epsilony);
		tauxy = E / (2 * (1 + v)) * gammaxy;
		Stress = sqrt(0.5 * ((sigmax - sigmay) * (sigmax - sigmay) + (sigmax - sigmaz) * (sigmax - sigmaz) + (sigmay - sigmaz) * (sigmay - sigmaz)) + 3 * tauxy * tauxy);
	}

	///////////////////////////////////////////
	//printing/drawing stuff
	void PrintStress() { std::cout << "element id: " << id << "stress: " << Stress << std::endl; }
	void PrintNodes()
	{
		std::cout << "wezly elementu id: " << id << " \n";
		for (auto it:nodes)
			std::cout << it;
	}
	void ColorElement(const double& min_stress, const double& max_stress)
	{
		T px[3] = { nodes[0].x, nodes[1].x, nodes[2].x };
		T py[3] = { nodes[0].y, nodes[1].y, nodes[2].y };
		setcolor((log(Stress) - min_stress) / (max_stress - min_stress)); //log
		polygon(px, py, 3);
		setgray(1);
	}

	///////////////////////////////////////////
	//getting access
	int GetDOF(const int& no_local_dof) const
	{
		int no_node = no_local_dof / 2;
		int dir = no_local_dof % 2;
		return nodes[no_node].id * 2 + dir;
	}
	T GetMass() const { return Mass; }
};

///////////////////////////////////////////
//class Mesh
template<class T>
class Mesh
{
private:
	std::vector<Node<T>> nodes;
	std::vector<Element<T>> elements;
	T max_element_dimension;
	int no_nodes = 0;
	int no_elements = 0;
	std::string mesh_name;
	template<class T> friend class ExactMassMatrix;
	template<class T> friend class ApproximateMassMatrix;
	template<class T> friend class StiffnessMatrix;
public:
	///////////////////////////////////////////
	//creating stuff
	Mesh(const T& max = 50) : max_element_dimension(max) {}
	~Mesh() = default;
	void AddNodes()
	{
		clear();
		double mar = 130, mag = 10, sx = 90, sy = 50;
		for (;;)
			if (mousedown())
				if (whichmousebutton() == 0)
				{
					T x = ((T)mouseclickx() - scale_dx - sx + mag - 1) / scale_x;
					T y = ((T)mouseclicky() - scale_dy - sy + mag - 1) / scale_y;
					Node<T> new_node(x, y, no_nodes);
					nodes.push_back(new_node);
					circle(new_node.x, new_node.y, 5);
					if (no_nodes > 0)
						line(nodes[no_nodes - 1].x, nodes[no_nodes - 1].y, nodes[no_nodes].x, nodes[no_nodes].y);
					no_nodes++;
				}
				else break;
		line(nodes[0].x, nodes[0].y, nodes[no_nodes - 1].x, nodes[no_nodes - 1].y);
	}
	void AddBoundaryCondition()
	{
		double mar = 130, mag = 10, sx = 90, sy = 50;
		T current_distance;
		int closest_node;
		for (;;)
			if (mousedown())
				if (whichmousebutton() == 0)
				{
					T x = ((T)mouseclickx() - scale_dx - sx + mag - 1) / scale_x;
					T y = ((T)mouseclicky() - scale_dy - sy + mag - 1) / scale_y;
					current_distance = nodes[0].CalculateDistance(x, y);
					closest_node = 0;
					for (auto it : nodes)
						if (it.CalculateDistance(x, y) < current_distance)
						{
							current_distance = it.CalculateDistance(x, y);
							closest_node = it.id;
						}
					nodes[closest_node].fix = true;
					setcolor(1);
					circle(nodes[closest_node].x, nodes[closest_node].y, 5);
					setgray(1);
				}
				else break;
	}
	void CreateMesh()
	{
		int i = 0, j = no_nodes - 1;
		for (;;)
		{
			elements.push_back(*new Element<T>(nodes[i], nodes[i + 1], nodes[j], no_elements));
			no_elements++;
			i++;
			if (no_elements == no_nodes - 2)
				break;
			elements.push_back(*new Element<T>(nodes[j], nodes[j - 1], nodes[i], no_elements));
			no_elements++;
			j--;
			if (no_elements == no_nodes - 2)
				break;
		}
		if (mesh_name == "car.txt")
			AddWheels();

		RefineMesh(max_element_dimension);
	}
	void AddWheels() //ONLY TO USE WITH CAR.TXT
	{
		nodes.push_back(*new Node<T>(62, 5, no_nodes));
		no_nodes++;
		nodes.push_back(*new Node<T>(81, 1, no_nodes, 0, 1));
		no_nodes++;
		nodes.push_back(*new Node<T>(100, 5, no_nodes));
		no_nodes++;

		nodes.push_back(*new Node<T>(326, 5, no_nodes));
		no_nodes++;
		nodes.push_back(*new Node<T>(345, 1, no_nodes, 0, 1));
		no_nodes++;
		nodes.push_back(*new Node<T>(364, 5, no_nodes));
		no_nodes++;

		elements.push_back(*new Element<T>(nodes[3], nodes[4], nodes[40], no_elements));
		no_elements++;
		elements.push_back(*new Element<T>(nodes[40], nodes[4], nodes[41], no_elements));
		no_elements++;
		elements.push_back(*new Element<T>(nodes[41], nodes[4], nodes[42], no_elements));
		no_elements++;
		elements.push_back(*new Element<T>(nodes[42], nodes[4], nodes[5], no_elements));
		no_elements++;

		elements.push_back(*new Element<T>(nodes[13], nodes[14], nodes[43], no_elements));
		no_elements++;
		elements.push_back(*new Element<T>(nodes[43], nodes[14], nodes[44], no_elements));
		no_elements++;
		elements.push_back(*new Element<T>(nodes[44], nodes[14], nodes[45], no_elements));
		no_elements++;
		elements.push_back(*new Element<T>(nodes[45], nodes[14], nodes[15], no_elements));
		no_elements++;
	}
	void RefineMesh(const T& max)
	{
		this->max_element_dimension = max;
		bool IsEveryElementOK = false;
		while (IsEveryElementOK == false)
		{
			IsEveryElementOK = true;
			DivideElements(IsEveryElementOK, max);
		}
		for (auto it = elements.begin(); it != elements.end(); it++)
			it->CalculateProperties();
	}
	void DivideElements(bool& IsEveryElementOK, const T& max)
	{
		const int no_old_elements = no_elements;
		bool ElementRefined = false;
		int node0 = 0, node1 = 1, node2 = 2;
		for (int i = 0; i < no_old_elements; i++)
		{
			if (elements[i].TheLongestSide(node0, node1, node2) > max)
			{
				IsEveryElementOK = false;
				ElementRefined = false;
				Node<T> new_node = elements[i].CalculateMidNode(node0, node1, no_nodes);
				for (auto it : nodes)
					if (it == new_node)
					{
						elements.push_back(*new Element<T>(it, elements[i].nodes[node1], elements[i].nodes[node2], no_elements));
						elements[i].nodes[node1] = it;
						no_elements++;
						ElementRefined = true;
					}
				if (ElementRefined == false)
				{
					nodes.push_back(new_node);
					elements.push_back(*new Element<T>(new_node, elements[i].nodes[node1], elements[i].nodes[node2], no_elements));
					elements[i].nodes[node1] = new_node;
					no_elements++;
					no_nodes++;
				}
			}
		}
	}

	///////////////////////////////////////////
	//calculating stuff
	void AddDisplacement(const std::vector<T>& vec)
	{
		for (int i = 0; i < no_nodes; i++)
		{
			nodes[i].x += vec[2 * i];
			nodes[i].y += vec[2 * i + 1];
		}
		for (int i = 0; i < no_elements; i++)
		{
			elements[i].nodes[0] = nodes[elements[i].nodes[0].id];
			elements[i].nodes[1] = nodes[elements[i].nodes[1].id];
			elements[i].nodes[2] = nodes[elements[i].nodes[2].id];
		}
	}
	void CalculateStresses(const std::vector<T>& vec)
	{
		for (auto it = elements.begin(); it != elements.end(); it++)
			it->CalculateStress(vec[it->nodes[0].id * 2], vec[it->nodes[0].id * 2 + 1], vec[it->nodes[1].id * 2], vec[it->nodes[1].id * 2 + 1], vec[it->nodes[2].id * 2], vec[it->nodes[2].id * 2 + 1]);
	}
	const T& MaxStress()
	{
		T max = 0;
		double scale = 0.8;
		for (auto it : elements)
			if (it.Stress > max)
				max = it.Stress;
		std::cout << "max stress: " << max << std::endl;
		return max * scale;
	}
	const T& MinStress()
	{
		T min = elements[0].Stress;
		for (auto it : elements)
			if (it.Stress < min)
				min = it.Stress;
		std::cout << "min stress: " << min << std::endl;
		return min;
	}

	///////////////////////////////////////////
	//printing/drawing stuff
	void DrawNodes()
	{
		clear();
		for (auto it:nodes)
			if (it.fix == 0)
				point(it.x, it.y);
			else
			{
				setcolor(1);
				circle(it.x, it.y, 2);
				setgray(1);
			}
	}
	void PrintNodes()
	{
		for (auto it:nodes)
			std::cout << it;
	}
	void ReadNodes(const char* file_name)
	{
		std::ifstream in;
		try
		{
			in.open(file_name);
			if (in.is_open() == false)
				throw std::string(file_name);
			mesh_name = file_name;
			in >> no_nodes;
			nodes.resize(no_nodes);
			int id;
			T x, y;
			bool fix;
			for (auto it = nodes.begin(); it != nodes.end(); it++)
			{
				in >> id >> x >> y >> fix;
				*it = *new Node<T>(x, y, id, fix);
			}
			in.close();
		}
		catch (std::string file_name)
		{
			std::cout << "file " << file_name << " not opened\n";
			exit(0);
		}
	}
	void WriteNodes(const char* file_name)
	{
		std::ofstream out;
		try
		{
			out.open(file_name);
			if (out.is_open() == false)
				throw std::string(file_name);
			out << no_nodes << "\n";
			for (auto it : nodes)
				out << it.id << "\t" << it.x << "\t" << it.y << "\t" << it.fix << std::endl;
			out.close();
		}
		catch (std::string file_name)
		{
			std::cout << "file " << file_name << " not opened\n";
			exit(0);
		}
	}
	void DrawMesh()
	{
		clear();
		DrawNodes();

		///SCENERY
		int road[8] = { 0, 375, 0, 500, 1200, 500, 1200, 375 };
		setfillstyle(1, LIGHTGRAY);
		fillpoly(4, road);

		int tree[8] = { 910, 375, 910, 200, 1010, 200, 1010, 375 };
		setfillstyle(1, BROWN);
		fillpoly(4, tree);

		setfillstyle(1, GREEN);
		fillellipse(960, 150, 100, 100);
		///

		for (auto it : elements)
		{
			line(it.nodes[0].x, it.nodes[0].y, it.nodes[1].x, it.nodes[1].y);
			line(it.nodes[1].x, it.nodes[1].y, it.nodes[2].x, it.nodes[2].y);
			line(it.nodes[2].x, it.nodes[2].y, it.nodes[0].x, it.nodes[0].y);
		}
	}
	void ColorMesh()
	{
		clear();

		///SCENERY
		int road[8] = { 0, 375, 0, 500, 1200, 500, 1200, 375 };
		setfillstyle(1, LIGHTGRAY);
		fillpoly(4, road);

		int tree[8] = { 910, 375, 910, 200, 1010, 200, 1010, 375 };
		setfillstyle(1, BROWN);
		fillpoly(4, tree);

		setfillstyle(1, GREEN);
		fillellipse(960, 150, 100, 100);
		///

		T min = log(MinStress());
		T max = log(MaxStress());
		for (auto it : elements)
			it.ColorElement(min, max);
	}
	void ColorMesh(const T& min, const T& max)
	{
		clear();

		///SCENERY
		int road[8] = { 0, 375, 0, 500, 1200, 500, 1200, 375 };
		setfillstyle(1, LIGHTGRAY);
		fillpoly(4, road);

		int tree[8] = { 910, 375, 910, 200, 1010, 200, 1010, 375 };
		setfillstyle(1, BROWN);
		fillpoly(4, tree);

		setfillstyle(1, GREEN);
		fillellipse(960, 150, 100, 100);
		///

		for (auto it : elements)
			it.ColorElement(log(5*min), log(0.2*max)); //5, 0.2 to color up a bit
	}
	void PrintStresses()
	{
		for (auto it : elements)
			it.PrintStress();
	}

	///////////////////////////////////////////
	//getting access
	const int& GetNoNodes() const { return no_nodes; }
	const int& GetNoElements() const { return no_elements; }
	const bool& GetBoundaryCondition(const int& id) const { return nodes[id].fix; }
	const bool& GetGroundCondition(const int& id) const { return nodes[id].ground; }
	int GetNodeX(const int& id) const { return nodes[id].x; }
	int GetNodeY(const int& id) const { return nodes[id].y; }
	const T& GetElementStress(const int& id) const { return elements[id].Stress; }
	void AddElementStress(const T& stress, const int& id) { elements[id].Stress = stress; }
};
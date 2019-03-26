#include "Trajectory_without_drag.h"
#include <iostream>
#include <iostream>
#include <fstream>

int main()
{
	Trajectory_without_drag trajectory;
	std::vector<Point3D> points;
	Point3D first, last;

	first.x = 1;
	last.x = 1.5;
	points.push_back(first);
	points.push_back(last);

	std::vector<Point3D> theta_ref = trajectory.calculate_theta_ref(points);



	std::ofstream myfile;
	myfile.open("theta_ref.csv");
	
	int number_of_points = theta_ref.size();
	std::cout << "Number of point is " << number_of_points << std::endl;
	for (auto& point : theta_ref) {

		
		
		myfile << point.x << "," << point.y << "," << point.z << ",";
		
		
		myfile << "\n";
		

	}
	
	myfile.close();
	std::cout << "Gotovo" << std::endl;

	std::cin.get();
}
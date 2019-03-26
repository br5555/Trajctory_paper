#pragma once
#include <vector>
#include <math.h>
#include <vector>

struct Point3D
{
	Point3D()
	{
		x = 0.0;
		y = 0.0;
		z = 0.0;
	}

	union
	{
		struct { double x, y, z; };
		struct { double theta_x, theta_y, theta_z; };
		
	};
	

	Point3D operator+ (Point3D const &obj) {
		Point3D res;
		res.x = x + obj.x;
		res.y = y + obj.y;
		res.z = z + obj.z;
		return res;
	}
};

constexpr int num_time_segments = 7;
class Trajectory_without_drag
{
public:
	Trajectory_without_drag();
	Trajectory_without_drag(double jerk_max, double  jerk_min, double  acc_max, double  acc_min, double vel_max, double  vel_min, double  vel_start, double  vel_end, double  p_start, double  p_end
		, double  acc_start, double  acc_end);	
	~Trajectory_without_drag();
	std::vector<Point3D> calculate_theta_ref(std::vector<Point3D> points);
	
private:
	bool calculate_trajectory_case1(double t[]);
	bool calculate_trajectory_case2(double t[]);
	bool calculate_trajectory_case3(double t[]);
	bool calculate_trajectory_case4(double t[]);
	bool calculate_trajectory_case5(double t[]);
	bool calculate_trajectory_case6(double t[]);
	bool calculate_trajectory_case7(double t[]);
	bool calculate_trajectory_case8(double t[]);
	inline double calculate_acceleration(double t, double acc_old, double jerk){
		return acc_old + t*jerk;
	}		
	inline double check_point()
	{	
		double tmp;
		if(m_p_start > m_p_end)
		{
			tmp = m_p_start;
			m_p_start = m_p_end;
			m_p_end = tmp;
			return -1.0;
		}
		else
			return 1.0;
	}	

	double m_jerk_max, m_jerk_min, m_acc_max, m_acc_min, m_vel_max, m_vel_min, m_vel_start, m_vel_end, m_p_start, m_p_end, m_acc_start, m_acc_end;
	void calculate_theta_for_two_points(char axis, std::vector<Point3D>& theta_refs);
};


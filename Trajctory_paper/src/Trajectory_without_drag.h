#pragma once
#include <vector>
#include <math.h>
#include <vector>

/// <summary>
/// Class representing point in 3D space
/// </summary>
struct Point3D
{
	/// <summary>
	/// Default constructor
	/// </summary>
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
	
	/// <summary>
	/// Overriding + operator
	/// </summary>
	/// <param name="obj">3D point</param>
	/// <returns>new 3D point</returns>
	Point3D operator+ (Point3D const &obj) {
		Point3D res;
		res.x = x + obj.x;
		res.y = y + obj.y;
		res.z = z + obj.z;
		return res;
	}
};
/// <summary>
/// Number of time segments
/// </summary>
constexpr int num_time_segments = 7;

/// <summary>
/// Class representing local planner from paper
/// author={M. {Beul} and S. {Behnke}}, 
/// booktitle = { 2016 International Conference on Unmanned Aircraft Systems(ICUAS) },
/// title = { Analytical time - optimal trajectory generation and control for multirotors },
/// </summary>
class Trajectory_without_drag
{
public:
	/// <summary>
	/// Default constructor
	/// </summary>
	Trajectory_without_drag();
	
	/// <summary>
	/// Constructor which defines local planner with physical limits
	/// </summary>
	/// <param name="jerk_max">maximum jerk</param>
	/// <param name="jerk_min">minimum jerk</param>
	/// <param name="acc_max">maximum acceleration</param>
	/// <param name="acc_min">minimum acceleration</param>
	/// <param name="vel_max">maximum velocity</param>
	/// <param name="vel_min">minimum velocity</param>
	/// <param name="vel_start">starting velocity</param>
	/// <param name="vel_end">ending velocity</param>
	/// <param name="p_start">starting point</param>
	/// <param name="p_end">ending point</param>
	/// <param name="acc_start">starting acceleration</param>
	/// <param name="acc_end">ending acceleration</param>
	Trajectory_without_drag(double jerk_max, double  jerk_min, double  acc_max, double  acc_min, double vel_max, double  vel_min, double  vel_start, double  vel_end, double  p_start, double  p_end, double  acc_start, double  acc_end);	
	
	/// <summary>
	/// Default destructor
	/// </summary>
	~Trajectory_without_drag();

	/// <summary>
	/// Calculating theta reference from desire set of points
	/// </summary>
	/// <param name="points"></param>
	/// <returns></returns>
	std::vector<Point3D> calculate_theta_ref(std::vector<Point3D> points);
	
private:
	/// <summary>
	/// Calculate trajectory using Case 1
	/// </summary>
	/// <param name="t">Calculated time intervals</param>
	/// <returns>True if trajectory is found otherwise false</returns>
	bool calculate_trajectory_case1(double t[]);

	/// <summary>
	/// Calculate trajectory using Case 2
	/// </summary>
	/// <param name="t">Calculated time intervals</param>
	/// <returns>True if trajectory is found otherwise false</returns>
	bool calculate_trajectory_case2(double t[]);

	/// <summary>
	/// Calculate trajectory using Case 3
	/// </summary>
	/// <param name="t">Calculated time intervals</param>
	/// <returns>True if trajectory is found otherwise false</returns>
	bool calculate_trajectory_case3(double t[]);

	/// <summary>
	/// Calculate trajectory using Case 4
	/// </summary>
	/// <param name="t">Calculated time intervals</param>
	/// <returns>True if trajectory is found otherwise false</returns>
	bool calculate_trajectory_case4(double t[]);

	/// <summary>
	/// Calculate trajectory using Case 5
	/// </summary>
	/// <param name="t">Calculated time intervals</param>
	/// <returns>True if trajectory is found otherwise false</returns>
	bool calculate_trajectory_case5(double t[]);

	/// <summary>
	/// Calculate trajectory using Case 6
	/// </summary>
	/// <param name="t">Calculated time intervals</param>
	/// <returns>True if trajectory is found otherwise false</returns>
	bool calculate_trajectory_case6(double t[]);

	/// <summary>
	/// Calculate trajectory using Case 7
	/// </summary>
	/// <param name="t">Calculated time intervals</param>
	/// <returns>True if trajectory is found otherwise false</returns>
	bool calculate_trajectory_case7(double t[]);

	/// <summary>
	/// Calculate trajectory using Case 8
	/// </summary>
	/// <param name="t">Calculated time intervals</param>
	/// <returns>True if trajectory is found otherwise false</returns>
	bool calculate_trajectory_case8(double t[]);

	/// <summary>
	/// Calculate discrete acceleration
	/// </summary>
	/// <param name="t">elapsed time</param>
	/// <param name="acc_old">accelration from previous discrete step </param>
	/// <param name="jerk">jerk at time t</param>
	/// <returns>acceleration at time t</returns>
	inline double calculate_acceleration(double t, double acc_old, double jerk){
		return acc_old + t*jerk;
	}

	/// <summary>
	/// Check points for sign (+ if next point is in positive direction from current point otherwise -)
	/// </summary>
	/// <returns>sign</returns>
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

	/// <summary>
	/// m_jerk_max- maximum jerk
	/// m_jerk_min - minimum jerk
	/// m_acc_max - maximum acceleration
	/// m_acc_min - minimum acceleration
	/// m_vel_max - maximum velocity
	/// m_vel_min - minimum velocity
	/// m_vel_start - starting velocity
	/// m_vel_end - ending velocity
	///  m_p_start - starting position 
	///  m_p_end - ending position 
	///  m_acc_start - starting acceleration
	///  m_acc_end - ending acceleration 
	/// </summary>
	double m_jerk_max, m_jerk_min, m_acc_max, m_acc_min, m_vel_max, m_vel_min, m_vel_start, m_vel_end, m_p_start, m_p_end, m_acc_start, m_acc_end;

	/// <summary>
	/// Calculate theta ref for sequential points
	/// </summary>
	/// <param name="axis"> x, y or z axis</param>
	/// <param name="theta_refs">theta reference </param>
	void calculate_theta_for_two_points(char axis, std::vector<Point3D>& theta_refs);
};


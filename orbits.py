from Vector import Vector3, Rotation
import sys
import copy
from math import *

def kepler1order(e, M, x):
    """First order kepler equation solver derived using taylor series
    expansion.
    e is orbit eccentricity.
    M is mean anomaly,
    x is previous Eccentric Anomaly"""
    return (x-e*sin(x)-M)/(1.0-e*cos(x))

def kepler2order(e, M, x):
    """Second order kepler equation solver derived using taylor series
    expansion.
    e is orbit eccentricity.
    M is mean anomaly,
    x is previous Eccentric Anomaly"""
    t1 = -1 + e*cos(x)
    t2 = e*sin(x)
    t3 = -x+t2+M
    return t3/(0.5*t3*t2/t1 + t1)

def solve_kepler(e, M):
    """An iterative function to solve kepler's equation."""
    E0 = M
    E = E0
    x0 = E0
    error = 1.0
    i = 0
    while(abs(error) > 1e-8):
        x1 = x0 - kepler1order(e, M, x0)
        E = x1
        error = x1 - x0
        x0 = x1

        #print("error: {0} E: {1}".format(error, E))

    return E

"""Both true anomaly and r vector make up the polar coordinates for
the space ship relative to the planet"""
def calc_true_anomaly_low_e(e, M):
    """This function is designed for orbits with low eccentricity, it is only
    accurate in those cases, but should be faster"""
    return M + 2.0*e*sin(M) + 1.25*(e**2) * sin(2.0*M)

def calc_true_anomaly(e, E):
    return acos((cos(E) - e)/(1.0-e*cos(E)))

def calc_r_vector(e, E, a):
    return a*(1.0 - e*cos(E))

def calc_orbit(vel, pos, planet):
    v = vel.mag()
    r = pos.mag()
    a = 1.0/(2.0/r - (v**2)/planet.u) #calculate semi major axis
    C = 2*planet.u/(r*(v**2))

def plot_vector_list2d(vectors, polygon=False, **kwargs):
    x = []
    y = []
    points = []
    for v in vectors:
        x.append(v.x)
        y.append(v.y)
        points.append([x, y])

    if polygon:
        pylab.Polygon(points, **kwargs)
    else:
        pylab.plot(x, y, **kwargs)

def plot_vector_arrow2d(origin, v, **kwargs):
    """plots a vector as an arrow. origin is the origin of the vector,
    and v is the vector"""
    head_length = 0.3 * v.mag()
    head_width = 0.08 * v.mag()

    end_pos = origin.add(v)
    v_norm = v.norm()
    end_left = v_norm.cross(Vector3(0,0,1)).fmul(head_width)
    end_right = end_left.neg()
    end_left = end_left.add(end_pos)
    end_right = end_right.add(end_pos)
    end_tip = end_pos.add(v_norm.fmul(head_length))

    plot_vector_list2d([origin, end_pos, end_left, end_tip, end_right, end_pos], **kwargs)



class PhysicsObject(object):
    G = 6.67e-11

    @classmethod
    def copy(cls, other):
        return cls(vel=Vector3.fromVector(other.vel),
            pos=Vector3.fromVector(other.pos),
            acc=Vector3.fromVector(other.acc),
            mass=other.mass,
            radius=other.radius,
            color=other.color,
            show_vel=other.show_vel)

    def __init__(self,**kwargs): #self,pos, vel, acc, mass, radius=None, color="grey", show_vel=True,):
        self.vel = kwargs['vel']
        self.pos = kwargs['pos']
        self.acc = kwargs['acc']
        self.mass = kwargs['mass']
        self.radius = kwargs['radius']
        self.color= kwargs['color']
        self.show_vel = kwargs['show_vel']

    def copy_motion(self, other):
        self.pos.copy(other.pos)
        self.vel.copy(other.vel)
        self.acc.copy(other.acc)

    def get_force(self):
        return Vector3()

    def get_torque(self):
        return Vector3()

    def calc_grav_constant(self):
        return self.mass*self.G
    
    def get_r(self, other):
        return self.pos.sub(other.pos)
    
    @classmethod
    def gravsum(cls, o1, o2):
        """Sum of the gravitational parameters of two bodies"""
        return cls.G * (o1.mass + o2.mass)
    
    @property
    def u(self):
        return self.calc_grav_constant()

    def plot(self, ax, threed=False):
        rad = 0.0
        if self.radius == None:
                if self.mass == 0:
                    rad = 1.0
                else:
                    rad = 1.0668117883456128e-18 * self.mass
        else:
            rad = self.radius

        if self.radius > 0:
            circle = pylab.Circle((self.pos.x,self.pos.y), color=self.color, radius=rad, clip_on=False)
            ax.add_patch(circle)

            #draw object crosshairs
            circle_vec_left = Vector3(-1.0, 0.0, 0.0).fmul(rad).add(self.pos)
            circle_vec_right = Vector3(1.0, 0.0, 0.0).fmul(rad).add(self.pos)
            circle_vec_up = Vector3(0.0, 1.0, 0.0).fmul(rad).add(self.pos)
            circle_vec_down = Vector3(0.0, -1.0, 0.0).fmul(rad).add(self.pos)
            plot_vector_list2d([circle_vec_left, circle_vec_right], color="red")
            plot_vector_list2d([circle_vec_down, circle_vec_up], color=self.color)

        plot_vector_list2d([self.pos], marker='o', color=self.color)

        #draw object velocity
        if(self.vel.mag() > 0.0 and self.show_vel):
            plot_vector_arrow2d(self.pos, self.vel.fmul(400), color="blue")
            """

            print(vel_end, self.pos)
            plot_vector_list2d([self.pos, vel_end])"""

    def __repr__(self):
        return "Physics Body [vel: {0} pos: {1} acc: {2} rad: {3} mass: {4}]".format(str(self.vel), str(self.pos), str(self.acc), str(self.radius), str(self.mass))


def angle_difference(targetA,sourceA):
    a = targetA - sourceA
    if a > 180:
        a -= 360
    if a < -180:
        a += 360
    return a;

def standard_angle_range(angle):
    pi2 = 2.0*pi
    if angle < 0.0:
        return angle + pi2

    if angle >= pi2:
        return angle - pi2

    return angle





class Orbit(object):
    def __init__(self, **args):
        if "color" in args.keys():
            self.color = args["color"]
        else:
            self.color =  "black"

class KeplerOrbit2(Orbit):
    @staticmethod
    def fromVectors(pos, vel, mu):
        """
        pos: position of element
        vel: velocity of element
        mu:  gravitational constant/parameter = GM
        """
        #calculate orbital momentum vector h (m^2/s)
        r = pos.mag()
        h = pos.cross(vel)
        #obtain eccentricity vector
        print("Vel Pos", vel, pos)
        e1 = vel.cross(h).fdiv(mu)
        e2 = pos.fdiv(pos.mag())
        print("e1 e2", e1, e2)
        e_vec = e1.sub(e2)
        e = e_vec.mag()
        print("evec", e_vec)
        #determine n (m^2/s)

        #calculate orbit inclination i by using orbital
        #momentum vector h, where h_z is third component of h
        i = acos(h.z/h.mag())
        while i > pi/2.0:
           i = i - pi

        print("i: ", i)


        n_vec = None
        if i == 0.0:
            n_vec = Vector3.X()
        else:
            n_vec = Vector3.Z().cross(h)

        #obtain longtitude of ascending node
        long_asc = 0.0
        if i == 0.0:
            long_asc = 0.0
        elif n_vec.y >= 0:
            print("n_vec:" , n_vec)
            print("h: ", h)
            long_asc = acos(n_vec.x/n_vec.mag())
        else:
            long_asc = 2.0*pi - acos(n_vec.x/n_vec.mag())

        #obtain argument of periapsis
        arg_per = 0.0
        if e_vec.z >= 0:
            arg_per = acos(n_vec.dot(e_vec)/(n_vec.mag() * e))
        else:
            arg_per = 2.0 * pi - acos(n_vec.dot(e_vec)/(n_vec.mag() * e))

        print("n_vec:" , n_vec)
        print("h: ", h)

        p = (h.mag()**2)/2.0 #semi-latus rectum

        a = 0.0
        semilatus = 0.0
        v = 0.0
        M = 0.0
        true_lon = 0.0
        arg_lat = 0.0
        lon_per = 0.0
        if e == 0:
            #circle
            a = r
            b = r

            lon_per = None
            arg_per = None

            if i==0.0:
                arg_lat = None
                if(vel.x > 0.0):
                    true_lon = 2*pi - acos(pos.x/r)
                else:
                    true_lon = acos(pos.x/r)
            else:
                if n_vec.dot(vel) > 0.0:
                    arg_lat = 2*pi - acos(n_vec.dot(pos)/(n_vec.mag()*r))
                else:
                    arg_lat = acos(n_vec.dot(pos)/(n_vec.mag()*r))

                #calculate true anomaly v (rad)
                if (pos.dot(vel) >= 0):
                    true_lon = Vector3(1,0,0).angle(pos)
                else:
                    true_lon = 2.0*pi - Vector3(1,0,0).angle(pos)

        elif (e > 0 and e < 1):
            #ellipse
            #calculate semi-major axis a
            a = 1.0/((2.0/pos.mag()) - ((vel.mag()**2)/mu))
            rp = (1 + e)*a #pericenter
            ra = (1 - e)*a #apocenter
            b = sqrt(ra*rp)  #semi-minor axis

            #calculate eccentric anomaly E
            E = 2 * atan(tan(v/2)/sqrt((1+e)/(1-e)))

            #calculate mean anomaly with the help of kepler's equation
            #from the eccentric anomaly, and the eccentricity of the orbit
            M = E - e * sin(E)

            #longtitude of periapsis
            long_per = long_asc + arg_per

            #calculate true anomaly v (rad)
            if (pos.dot(vel) >= 0):
                v = e_vec.angle(pos)
            else:
                v = 2.0*pi - e_vec.angle(pos)

            #argument of latitude
            arg_lat = v + arg_per
            #true longtitude
            true_lon = v + long_per

        elif e == 1:
            #parabola
            b = sys.float_info.max
            a = sys.float_info.max #infinity

        elif e > 1:
            #hyperbola
            a = -p/(e**2 - 1.0)
            b = p/sqrt(e**2 - 1)



        print("True anomaly: ", degrees(v))

        return KeplerOrbit2(a,
                b,
                e,
                i,
                long_asc,
                arg_per,
                M,
                true_lon,
                arg_lat,
                lon_per,
                p,
                mu,
                n_vec,
                h)

    @staticmethod
    def calc_mean_motion(T):
        #return sqrt(6.67384e-11 * (self.parentbody.mass + self.childbody.mass) / (self.a**2))
        #return sqrt(6.67e-11 * (self.parentbody.mass) / (self.a**2))
        return 2*pi/T

    def kepler1order(self, M, x):
        """First order kepler equation solver derived using taylor series
        expansion.
        e is orbit eccentricity.
        M is mean anomaly,
        x is previous Eccentric Anomaly
        returns E eccentric anomaly"""
        return (x- self.e * sin(x)-M)/(1.0 - self.e * cos(x))

    def kepler2order(self, M, x):
        """Second order kepler equation solver derived using taylor series
        expansion.
        e is orbit eccentricity.
        M is mean anomaly,
        x is previous Eccentric Anomaly
        returns E eccentric anomaly"""
        t1 = -1 + self.e * cos(x)
        t2 = self.e * sin(x)
        t3 = -x+t2+M
        return t3/(0.5*t3*t2/t1 + t1)

    def solve_kepler(self, M):
        """An iterative function to solve kepler's equation."""
        E0 = M
        E = E0
        x0 = E0
        error = 1.0
        i = 0
        while(abs(error) > 1e-8):
            x1 = x0 - self.kepler2order(M, x0)
            E = x1
            error = x1 - x0
            x0 = x1

            #print("error: {0} E: {1}".format(error, E))
            i += 1

        return E #this is a hack, kepler solver isn't working properly!!

    def calc_mean_anomaly(self, dt):
        print("n:", self.n)
        return self.n * dt
        #E = self.calc_eccentric_anomaly()
        #M0 = E - self.e * sin(E)
        #return M0 + self.n*dt

    def calc_true_anomaly(self, E):
        #v = acos((cos(E) - self.e)/(1.0-e*cos(E)))
        #if E > pi:
        #    v = 2.0 * pi - v
        #v = 2.0 * atan(sqrt((1.0 + self.e)/(1.0 - self.e)) * tan(E/2.0))
        #if v < 0.0:
        #    v += 2.0*pi
        v = 2.0*atan2(sqrt(1+self.e)*sin(E/2.0),sqrt(1-self.e)*cos(E/2.0))
        return v

    def calc_r_from_true_anomaly(self, true_anomaly):
        #return self.p / (1.0 + self.e * cos(true_anomaly))
        return self.a * (1.0 - self.e**2)/(1.0 + self.e * cos(true_anomaly))

    def calc_r_from_E(self, E):
        return self.a * (1.0 - self.e * cos(E))

    def calc_pos_from_E(self, E):
        #todo: can optimise this a lot
        r = self.calc_r_from_E(E)
        v = self.calc_true_anomaly(E)
        angle_sum = self.arg_per + v

        x = Vector3.X()
        n_dir = x.rotate(Rotation.aroundZ(self.long_asc))
        h_dir = self.h_vec.norm()
        r_dir = n_dir.rotate(Rotation.aroundVector(h_dir, v))
        r_vec = r_dir.fmul(self.calc_r_from_E(E))
        r_vec = r_vec.rotate(Rotation.aroundVector(h_dir, self.arg_per))
        return r_vec


        #x = r * (cos(self.long_asc) * cos(angle_sum) - sin(self.long_asc) * sin(angle_sum) * cos(self.i))
        #y = r * (sin(self.long_asc) * cos(angle_sum) + cos(self.long_asc) * sin(angle_sum) * cos(self.i))
        #z = r * (sin(self.i)*sin(angle_sum))
        #return Vector3(x, y, z)

    def to_vectors(self):
        return

    def __init__(self,
                a,
                b,
                e,
                i,
                long_asc,
                arg_per,
                M,
                true_lon,
                arg_lat,
                lon_per,
                p,
                mu,
                n_vec,
                h_vec,
                **args):

        """
        a:         semi major axis (km)
        b:         semi minor axis (km)
        e:     eccentricity (unitless)
        i:         inclination (radians)
        long_asc: longtitude of ascending node (radians)
        arg_per:   argument of perigee (radians)
        M:         mean anomaly (radians)
        true_lon:  True Longtitude (radians)
        arg_lat:   argument of latitude (radians)
        long_per:  Longtitude of periapse (radians)
        p: semilatus rectum (km)
        mu:  gravitational constant
        n: ???
        """
        super(KeplerOrbit2, self).__init__(**args)
        self.a = a
        self.b = b
        self.e = e
        self.i = i
        self.long_asc = long_asc
        self.arg_per = arg_per
        self.M = M
        self.true_lon = true_lon
        self.arg_lat = arg_lat
        self.lon_per = lon_per
        self.p = p
        self.mu = mu
        self.n_vec = n_vec
        self.h_vec = h_vec

        if e == 0:
            print("Circular")
            self.c = 0.0
            self.ra = a
            self.rb = b
        elif (e > 0 and e < 1):
            #ellipse
            print("Ellipse")
            self.rp = (1 + e)*a #pericenter
            self.ra = (1 - e)*a #apocenter
            self.c = a * e #focal distance
        elif e == 1:
            #parabola
            print("Parabola")
            self.rp = p/2.0 #pericenter
            self.ra = sys.float_info.max
            self.c = sys.float_info.max
        elif e > 1:
            #hyperbola
            print("Hyperbola")
            self.rp = p/(1.0 + e)
            self.ra = None
            self.c = a * e

        #calculate mean motion
        #self.n = self.calc_mean_motion(self.T)

        print("a: ", self.a, "b: ", self.b)
        print("long_asc: ", degrees(self.long_asc))
        print("Argument of Periapsis: ", degrees(self.arg_per))
        print("rp: ", self.rp, "ra: ", self.ra)
        print("e: ", self.e)
        print("i: ", self.i)

        #calculate orbit period
        self.T = 2.0*pi*sqrt((a**3)/mu)
        self.n_vec = n_vec.fmul(self.calc_mean_motion(self.T))
        self.n = self.n_vec.mag()

        print("T: ", self.T)
        print("n: ", self.n)

    @staticmethod
    def fromOrbitalElements(a, b, e, i, long_asc, arg_per, M, true_lon, arg_lat, lon_per, p, mu):
        """
        a:         semi major axis (km)
        b:         semi minor axis (km)
        e:     eccentricity (unitless)
        i:         inclination (radians)
        long_asc: Right ascention of the ascending node (radians)
        arg_per:   argument of perigee (radians)
        M:         mean anomaly (radians)
        true_lon:  True Longtitude (radians)
        arg_lat:   argument of latitude (radians)
        long_per:  Longtitude of periapse (radians)
        p: semilatus rectum (km)
        mu:  gravitational constant
        """

        hval = b*sqrt(mu/a)
        h_vec1 = Vector3.Z().rotate(Rotation.aroundY(i))
        h_vec = h_vec1.rotate(Rotation.aroundZ(long_asc - pi/2.0))
        x = Vector3.X()
        rot_long_asc = Rotation.aroundZ(long_asc)
        n_dir = x.rotate(rot_long_asc).norm()
        n_vec = Vector3.Z().cross(h_vec)

        KeplerOrbit2(
                a,
                b,
                e,
                i,
                long_asc,
                arg_per,
                M,
                true_lon,
                arg_lat,
                lon_per,
                p,
                mu,
                n_vec,
                h_vec)

    def solve_dt(self, dt):
        M = self.calc_mean_anomaly(dt)
        M = standard_angle_range(M)
        E = self.solve_kepler(M)
        v = self.calc_true_anomaly(E)

        return [M, E, v]

    def get_orbit_period(self):
        return self.T

    def get_flight_angle(self):
        return self.flight_angle

    def get_periapse(self):
        return self.rp

    def get_apoapse(self):
        return self.ra


    def plot(self, ax):
        """
        """
        pos_list = []
        step_size = self.T/100 #step size proportial to the total orbit period.
        #need to make step size proportional to maximum velocity as well
        #multiply time by velocity to get distance based step size
        #sort of done that, variable names could be better
        print("step_size: ", step_size)
        dt = 0.0
        while dt < self.T:
            M = self.calc_mean_anomaly(dt)
            M = standard_angle_range(M)
            E = self.solve_kepler(M)
            #pos_old = self.calc_pos_old(E) #todo: use positions to plot orbit
            pos = self.calc_pos_from_E(E)
            radius = pos.mag()
            step_ratio = radius/self.get_apoapse()
            print(pos)

            print("step_ratio", step_ratio)
            print("Mean Anomaly: ", degrees(M))
            print("Eccentric Anomaly: ", degrees(E))
            print("True Anomaly: ", degrees(self.calc_true_anomaly(E)))
            print("Anomaly Diff: ", angle_difference(degrees(self.calc_true_anomaly(E)), degrees(E)))
            #print("Radius method 1: ", self.calc_r_from_true_anomaly(self.calc_true_anomaly(E)))
            #print("Radius method 2: ", self.calc_r_from_E(E))
            #print("Radius method 3: ", pos.mag())
            #print("Velocity: ", self.calc_vel_from_E(E))
            print("")

            pos_list.append(pos)
            #circle = pylab.Circle((pos.x, pos.y), radius=4e5, color=(sc, sc, sc), clip_on=False)
            #ax.add_patch(circle)
            dt += step_size * step_ratio
            #break


        #pos_list.append(pos_list[0]) # link back to start
        #print(pos_list[5])
        plot_vector_list2d(pos_list, color="orange")

class ProgressIndicator(object):
    def __init__(self, start_val, end_val):
        self.start_val = start_val
        self.end_val = end_val
        self.percentage_old = 0.0
        self.percentage_new = 0.0
        
    def update(self, new_val):
        if(new_val > 0.0):
            self.percentage_new = int((new_val/self.end_val)*100)
            
        if(self.percentage_new != self.percentage_old):
            print("Calculating: {0}%".format(self.percentage_new))
            self.percentage_old = self.percentage_new

class NumericOrbit(Orbit):
    def __init__(self, m_big, m_small, **args):
        super(NumericOrbit, self).__init__(**args)
        self.m_big = m_big
        self.m_small = m_small
        
    def get_energy(self, o1, o2):
        """gets the energy of o2 orbiting around o1"""
        u = PhysicsObject.gravsum(o1, o2) 
        e = o2.vel.mag()**2 - u/(o2.get_r(o1).mag())
        return e
    
    def energy_vel_adj(self, o1, o2, target_energy):
        """Adjusts the velocity of o2 to reach a target
        orbital energy orbiting around o1"""
        
        curr_energy = self.get_energy(o1, o2)
        pass #need to implement a root finding algorithm,
        #probably would be slower than just using more samples!
        
    def f_o1_on_o2(self, o1, o2):
        r_vec = o2.pos.sub(o1.pos)
        rhat_vec = r_vec.norm()
        dist = r_vec.mag()
        f_mag = -PhysicsObject.G * (o1.mass * o2.mass)/(dist**2)
        return rhat_vec.fmul(f_mag)

    def get_force(self, status):
        force = Vector3()
        force = force.add(self.f_o1_on_o2(self.m_big, status))
        force = force.add(status.get_force())
        return force
    
    def acceleration_o2_o1(self, o1, o2):
        """The "acceleration" method computes the acceleration between two
         bodies. It does not include mass in the calculation for simplicty.
         it is broken"""
        r_vec = o2.pos.sub(o1.pos)
        rhat_vec = r_vec.norm()
        dist = r_vec.mag()
        a_mag = -PhysicsObject.G * o1.mass/(dist**2)
        return rhat_vec.fmul(a_mag)
    
    def acceleration_mass(self, o1, o2):
        return self.f_o1_on_o2(o1, o2).fdiv(o2.mass)
    
    def calc_orbit(self):
        return []
    
    def plot(self, ax):
        pos_list = []
        t_list = self.calc_orbit(100000, 100)
        
        for status in t_list:
            pos_list.append(status.pos)
        
        plot_vector_list2d(pos_list, color=self.color)


class EulerOrbit(NumericOrbit):
    def __init__(self, m_big, m_small, **args):
        super(EulerOrbit, self).__init__(m_big, m_small, **args)
        
        if "conserve_energy" in args.keys():
            self.conserve_energy = args["conserve_energy"]
        else:
            self.conserve_energy = False
    
    def calc_orbit(self, t, t_step):
        t_list = []
        t_list.append(copy.deepcopy(self.m_small))
        
        last_status = t_list[0]
        
        t_curr = 0.0
        print("Calculating Euler Orbit")
        progress = ProgressIndicator(0.0, t)
        while t_curr < t:
            #for multibody gravity, just sum up all the 
            #forces right here
            new_f = self.f_o1_on_o2(self.m_big, last_status)
            new_a = new_f.fdiv(last_status.mass)
            new_status = copy.deepcopy(last_status)
            new_status.vel = new_status.vel.add(new_a.fmul(t_step))
            new_status.pos = new_status.pos.add(new_status.vel.fmul(t_step))
            t_list.append(new_status)
            progress.update(t_curr)
            last_status = new_status
            t_curr += t_step
            
        return t_list
        

class RK4Orbit(NumericOrbit):
    def __init__(self, m_big, m_small, **args):
        super(RK4Orbit, self).__init__(m_big, m_small, **args)
        
    def derivative(self, o):
        pass
    
    def calc_orbit(self, t, t_step):
        t_list = []
        t_list.append(copy.deepcopy(self.m_small))
        
        last_status = t_list[0]
        
        t_curr = 0.0
        #print("Calculating Euler Orbit")
        #progress = ProgressIndicator(0.0, t)
        i = 0
        while t_curr < t:
            k1_status = last_status
            k1_status.acc =  self.f_o1_on_o2(self.m_big, k1_status).fdiv(k1_status.mass)
            
            k2_status = copy.deepcopy(last_status)
            k2_status.pos = last_status.pos.add(last_status.vel.fmul(t_step*0.5))
            k2_status.vel = last_status.vel.add(k1_status.acc.fmul(t_step*0.5))
            k2_status.acc = self.f_o1_on_o2(self.m_big, k2_status).fdiv(k2_status.mass)
            
            k3_status = copy.deepcopy(last_status)
            k3_status.pos = last_status.pos.add(k2_status.vel.fmul(t_step*0.5))
            k3_status.vel = last_status.vel.add(k2_status.acc.fmul(t_step*0.5))
            k3_status.acc = self.f_o1_on_o2(self.m_big, k3_status).fdiv(k3_status.mass)
            
            k4_status = copy.deepcopy(last_status)
            k4_status.pos = last_status.pos.add(k3_status.vel.fmul(t_step))
            k4_status.vel = last_status.vel.add(k3_status.acc.fmul(t_step))
            k4_status.acc = self.get_force(k4_status).fdiv(k4_status.mass)
            
            
            new_status = copy.deepcopy(last_status)
            # xf = x + (dt/6.0)*(v1 + 2*v2 + 2*v3 + v4)
            new_status.pos = last_status.pos.add(k1_status.vel.add(k2_status.vel.fmul(2.0)).add(k3_status.vel.fmul(2.0)).add(k4_status.vel).fmul(t_step/6.0))
            # vf = v + (dt/6.0)*(a1 + 2*a2 + 2*a3 + a4)
            new_status.vel = last_status.vel.add(k1_status.acc.add(k2_status.acc.fmul(2.0)).add(k3_status.acc.fmul(2.0)).add(k4_status.acc).fmul(t_step/6.0))
            
            t_list.append(new_status)
            #progress.update(t_curr)
            last_status = new_status
            
            i += 1
            t_curr = t_step * i #instead of adding, to increase accuracy
            
            
        return t_list
        
        
class VelocityVerletOrbit(EulerOrbit):
    def __init__(self, m_big, m_small, **args):
        super(VelocityVerletOrbit, self).__init__(m_big, m_small, **args)
    
    def calc_orbit(self, t, t_step):
        t_list = []
        t_list.append(copy.deepcopy(self.m_small))
        
        last_status = t_list[0]
        last_status.acc = self.acceleration_o2_o1(self.m_big, last_status)
        
        t_curr = 0.0
        print("Calculating Verlet Orbit")
        progress = ProgressIndicator(0.0, t)
        i = 0
        while t_curr < t:
            new_status = copy.deepcopy(last_status)
            
            new_status.pos = new_status.pos.add(last_status.vel.fmul(t_step)).add(last_status.acc.fmul(0.5 * t_step**2))
            new_status.acc = self.acceleration_o2_o1(self.m_big, new_status)
            asum = last_status.acc.add(new_status.acc)
            new_status.vel = new_status.vel.add(asum.fmul(0.5 * t_step))
            
            t_list.append(new_status)
            progress.update(t_curr)
            last_status = new_status
            
            i += 1
            t_curr = t_step * i #instead of adding, to increase accuracy
            
            
        return t_list
    
class VerletOrbit(EulerOrbit):
    """Time Corrected Verlet Orbit, varying timestep"""
    def __init__(self, m_big, m_small, **args):
        super(VerletOrbit, self).__init__(m_big, m_small, **args)
    
    def calc_orbit(self, t, t_step):
        t_list = []
        t_list.append(copy.deepcopy(self.m_small))
        
        last_status = t_list[0]
        last_status.acc = self.acceleration_o2_o1(self.m_big, last_status)
        
        
        t_curr = 0.0
        print("Calculating Verlet Orbit")
        progress = ProgressIndicator(0.0, t)
        i = 0
        while t_curr < t:
            new_status = copy.deepcopy(last_status)
            
            if(i == 0):
                new_status.pos = new_status.pos.add(last_status.vel.fmul(t_step)).add(last_status.acc.fmul(0.5 * t_step**2))
            else:
                new_status.pos = t_list[i].pos.fmul(2.0).sub(t_list[i-1].pos).add(t_list[i].acc.fmul(t_step**2))
                
            new_status.acc = self.acceleration_o2_o1(self.m_big, new_status)
            asum = last_status.acc.add(new_status.acc)
            new_status.vel = new_status.vel.add(asum.fmul(0.5 * t_step))
            
            t_list.append(new_status)
            progress.update(t_curr)
            last_status = new_status
            
            i += 1
            t_curr = t_step * i #instead of adding, to increase accuracy
            
        return t_list
    
class TCVerletOrbit(EulerOrbit):
    """Time Corrected Verlet Orbit, varying timestep"""
    def __init__(self, m_big, m_small, **args):
        super(TCVerletOrbit, self).__init__(m_big, m_small, **args)
    
    def calc_orbit(self, t, t_step):
        t_list = []
        t_list.append(copy.deepcopy(self.m_small))
        
        last_status = t_list[0]
        last_status.acc = self.acceleration_o2_o1(self.m_big, last_status)
        
        t_step_adj = 2
        
        t_curr = 0.0
        t_step_curr = (t_step_adj * t_step)/last_status.acc.mag()
        t_step_last = t_step_curr
        print("Calculating Time Corrected Verlet Orbit")
        progress = ProgressIndicator(0.0, t)
        i = 0
        while t_curr < t:
            new_status = copy.deepcopy(last_status)
            
            if(i == 0):
                new_status.pos = new_status.pos.add(last_status.vel.fmul(t_step_curr)).add(last_status.acc.fmul(0.5 * t_step_curr**2))
            else:
                new_status.pos = t_list[i].pos.add(t_list[i].pos.sub(t_list[i-1].pos).fmul(t_step_curr/t_step_last)).add(t_list[i].acc.fmul(t_step_curr**2))
                
            new_status.acc = self.acceleration_o2_o1(self.m_big, new_status)
            asum = last_status.acc.add(new_status.acc)
            new_status.vel = new_status.vel.add(asum.fmul(0.5 * t_step))
            
            t_list.append(new_status)
            progress.update(t_curr)
            last_status = new_status
            
            t_step_last = t_step_curr
            t_step_curr = (t_step_adj * t_step)/new_status.acc.mag()
            
            i += 1
            t_curr += t_step_curr
            
        return t_list
        
        

def test_kepler_orbit():
    ##Creates empty axes (aspect=1 means scale things so that circles look like circles)
    ax = pylab.axes(aspect=1)
    earth = PhysicsObject(Vector3(), Vector3(), Vector3(),  5.972E24, 6371e3, show_vel = False)
    earth.plot(ax)

    iss2 = PhysicsObject(Vector3(0, float(6371e3 + 418e3), 0),
                        Vector3(-7.66e3*1.0, -1e3*3, 0), Vector3(),
                        450000.0, 0.0, color="green")

    iss = PhysicsObject(Vector3(float(6371e3 + 418e3), 0, 0),
                        Vector3(0, 7.66e3*1.2, 0), Vector3(),
                        450000.0, 0.0, color="green")
#
#     iss = PhysicsObject(Vector3(0.0, float(6371e3 + 418e3), 0),
#                         Vector3(7.66e3*1.2, 0, 0),
#                         0.0, 0.0, color="green")

    #iss_orbit1 = KeplerOrbit(earth, iss)
    #iss_orbit1.plot(ax)
    #print("\n\nORBIT 2:\n")
    #iss_orbit2 = KeplerOrbit2.fromVectors(iss.pos, iss.vel, earth.u)
    #iss_orbit2.plot(ax)
    
    
    #iss_orbit4 = VerletOrbit(earth, iss2)
    #iss_orbit4.plot(ax)
    
    #iss_orbit4 = TCVerletOrbit(earth, iss2, color="green")
    #iss_orbit4.plot(ax)
    
    iss_orbit4 = RK4Orbit(earth, iss2)
    iss_orbit4.plot(ax)
    
    print("\n\nORBIT 3:\n")
    iss_orbit3 = KeplerOrbit2.fromVectors(iss2.pos, iss2.vel, earth.u)
    iss_orbit3.plot(ax)
    iss2.plot(ax)

    pylab.show()

def test_kepler_solver():
    """expected output: True anomaly of around 151.28"""
    E = solve_kepler(0.1, 2.53755)
    print("True Anomaly: ", degrees(calc_true_anomaly(0.1, E)))

"""We want to find a specific e and M value given a position and
a velocity relative to the planet, plus the mass of planet and mass of object"""


"""In the future it would be good to scale everything to fit within a 1x1x1 box for the numerical
methods to have better accuracy. rescale to display them"""

if __name__ == "__main__":
    test_kepler_orbit()

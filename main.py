from kivy.app import App
from kivy.uix.widget import Widget
from kivy.uix.label import Label
from kivy.graphics import Color, Ellipse, Line, Triangle, Rectangle, Point
from kivy.clock import Clock
from kivy.config import Config
from orbits import PhysicsObject, RK4Orbit
from Vector import Vector3
import time
import threading
import copy
import random
import sys

if sys.version_info.major == 2:
	from Queue import Queue
else:
	from queue import Queue

class Planet(PhysicsObject):
	def __init__(self, **kwargs):
		super(Planet, self).__init__(**kwargs)

class Ship(PhysicsObject):
	def __init__(self, **kwargs):
		super(Ship, self).__init__(**kwargs)
		
		if 'direction' in kwargs:
			self.direction = kwargs['direction']
		else:
			self.direction = Vector3.Y()

		self.engine_power = 0.0
		self.engine_rating = 1e6
		self.isp = 1e6
		self.fuel = 100.0

	def copy_motion(self, other):
		super(Ship, self).copy_motion(other)
		self.engine_power = other.engine_power
		self.direction = other.direction

	def get_effective_engine_power(self):
		if self.fuel > 0.0:
			return self.engine_power
		else:
			return 0.0

	def get_engine_force_mag(self):
		return self.get_effective_engine_power() * self.engine_rating

	def get_engine_force(self):
		return self.direction.fmul(self.get_engine_force_mag())

	def get_force(self):
		return self.get_engine_force()

	def get_fuel_usage(self):
		return self.get_engine_force_mag() / self.isp

	def update(self, dt):
		fuel_dt = self.get_fuel_usage() * dt
		if fuel_dt > 0.0:
			if self.fuel - fuel_dt >= 0:
				self.fuel -= fuel_dt
			else:
				self.fuel = 0.0
			print(self.fuel)

class OrbitPath(Queue):
	def __init__(self, physics_object):
		Queue.__init__(self)
		self.physics_object = physics_object
		self.mutex = threading.Lock()

	def reset(self):
		while len(self) > 0:
			self.get()

	def reset_back(self, last_i):
		self.mutex.acquire()
		self.queue[(last_i+1):] = []
		self.mutex.release()

	def __len__(self):
		return self.qsize()

	def __iter__(self):
		return self.queue.__iter__()

	def put_items(self, t, dt, t_obj):
		self.put((t,dt,t_obj))

	def put(self, item):
		self.mutex.acquire()
		Queue.put(self, item)
		self.mutex.release()

	def peek(self, i):
		return self.queue[i]

	def join_t_list(self, new_t_list):
		self.t_list += new_t_list

	def get_min_t(self):
		min_t = sys.float_info.max
		min_t_status = None
		self.mutex.acquire()
		for status in self:
			if status[0] < min_t:
				min_t_status = status
				min_t = status[0]

		self.mutex.release()
		return min_t_status

	def get(self, blocking=True, timout=True):
		self.mutex.acquire()
		Queue.get(self, blocking, timout)
		self.mutex.release()

	def get_n(self, n, blocking=True, timeout=None):
		l = []
		i = 0
		while i < n:
			l.append(self.get(blocking, timeout))
			i+=1

		return l

	def get_status_interp(self, t, trim=True):
		min_t_status = self.get_min_t()
		if t <= min_t_status[0]:
			return min_t_status[2]

		i = 1
		prev_status = self.peek(0)
		while i < len(self):
			next_status = self.peek(i)

			if self.t_between_status(t, prev_status, next_status):
				retval =  self.interp_between_status(t, prev_status, next_status)
				#remove the first item from the queue, closest to the body
				if trim and (i > 1):
					self.get_n(i)

				return retval

			prev_status = next_status
			i+=1

	def t_between_status(self, t, status1, status2):
		if t > status1[0]:
			if t < status2[0]:
				return True

		return False

	def interp_between_status(self, t, status1, status2):
		retval = Ship.copy(self.physics_object)
		obj1 = status1[2]
		obj2 = status2[2]
		dt = status2[0] - status1[0]
		t_from_1 = t - status1[0]
		t_frac = t_from_1/dt
		dp = obj2.pos.sub(obj1.pos)
		dv = obj2.vel.sub(obj1.vel)

		p_from_1 = dp.fmul(t_frac) #linearly interpolated position from status1
		interp_p = obj1.pos.add(p_from_1)

		v_from_1 = dv.fmul(t_frac) #linearly interpolated velocity from status1
		interp_v = obj1.vel.add(v_from_1)

		retval.pos = interp_p
		retval.vel = interp_v
		retval.acc = obj1.acc

		return retval

	def get_point_list(self, center_vec, world_scale, step=1):
		point_list = []
		i = 0
		while i < (len(self)-(step-1)):
			status = self.peek(i)
			t = status[0]
			dt = status[1]
			phys_obj = status[2]
			pos = center_vec.add(phys_obj.pos.fmul(world_scale))
			point_list.append(pos.x)
			point_list.append(pos.y)
			i+=step

		return point_list

	def __repr__(self):
		self.mutex.acquire()
		s = "OrbitPath: ["
		for t, dt, obj in self.queue:
			s += "(t:{0},dt:{1},PhysicsObject[pos: {2}, vel: {3}]),\n".format(t,dt,obj.pos,obj.vel)

		s+="]"
		self.mutex.release()
		return s

class PreviewOrbitThread(threading.Thread):
	def __init__(self, t_max, dt, orbit_path, orbit, start_i, preview):
		super(PreviewOrbitThread, self).__init__()
		self.orbit = orbit
		self.orbit_path = orbit_path
		self.start_i = start_i
		self.t_max = t_max
		self.dt = dt
		self.preview = preview

	def stop(self):
		self.orbit.stop()

	def stopped(self):
		return self.orbit.stopped()

	def run(self):
		self.orbit.calc_orbit(self.t_max, self.dt, self.orbit_path, self.start_i, self.preview)

class RK4PreviewOrbit(RK4Orbit):
	def __init__(self, m_big, m_small, **args):
		super(RK4PreviewOrbit, self).__init__(m_big, m_small, **args)
		self._stop = threading.Event()
		self.thread = None

		self.k1_status = copy.deepcopy(m_small)
		self.k2_status = copy.deepcopy(m_small)
		self.k3_status = copy.deepcopy(m_small)
		self.k4_status = copy.deepcopy(m_small)

	def stop(self):
		self._stop.set()

	def stopped(self):
		return self._stop.isSet()

	def clear_stopped(self):
		self._stop.clear()

	def get_force_preview(self, status):
		return self.f_o1_on_o2(self.m_big, status)

	def get_force_realtime(self, status):
		force = Vector3()
		force = force.add(self.f_o1_on_o2(self.m_big, status))
		force = force.add(self.m_small.get_force())
		return force


	def get_force(self, status, preview, i):
		if i == 0 and not preview:
			return self.get_force_realtime(status)
		else:
			return self.get_force_preview(status)

	def stop_wait(self):
		if self.thread is not None:
			if not self.thread.stopped():
				self.thread.stop()
				self.thread.join()
				self.clear_stopped()

	def triggger_calc_orbit(self, t_max, dt, orbit_path, reset=True, preview=False):
		start_i = 0
		if reset:
			start_i = 0
			self.stop_wait()
			orbit_path.reset()
			reset_obj = copy.deepcopy(orbit_path.physics_object)
			orbit_path.put_items(0.0, dt, reset_obj)
		else:
			start_i = len(orbit_path) - 1

		self.stop_wait()
		self.thread = PreviewOrbitThread(t_max, dt, orbit_path, self, start_i, preview)

		self.thread.daemon = True
		self.thread.start()


	def calc_orbit(self, t_max, t_step, orbit_path, start_i, preview):
		if len(orbit_path) > 1:
			orbit_path.reset_back(start_i)

		last_status = orbit_path.peek(start_i)
		last_status_obj = last_status[2]

		t_curr = last_status[0]
		i = 0
		while t_curr < t_max:
			time.sleep(0.01)
			if self.stopped():
				break

			self.k1_status.copy_motion(last_status_obj)
			self.k1_status.acc =  self.f_o1_on_o2(self.m_big, self.k1_status).fdiv(self.k1_status.mass)

			self.k2_status.copy_motion(last_status_obj)
			self.k2_status.pos = last_status_obj.pos.add(last_status_obj.vel.fmul(t_step*0.5))
			self.k2_status.vel = last_status_obj.vel.add(self.k1_status.acc.fmul(t_step*0.5))
			self.k2_status.acc = self.f_o1_on_o2(self.m_big, self.k2_status).fdiv(self.k2_status.mass)

			self.k3_status.copy_motion(last_status_obj)
			self.k3_status.pos = last_status_obj.pos.add(self.k2_status.vel.fmul(t_step*0.5))
			self.k3_status.vel = last_status_obj.vel.add(self.k2_status.acc.fmul(t_step*0.5))
			self.k3_status.acc = self.f_o1_on_o2(self.m_big, self.k3_status).fdiv(self.k3_status.mass)

			self.k4_status.copy_motion(last_status_obj)
			self.k4_status.pos = last_status_obj.pos.add(self.k3_status.vel.fmul(t_step))
			self.k4_status.vel = last_status_obj.vel.add(self.k3_status.acc.fmul(t_step))
			k4_force = self.get_force(self.k4_status, preview, i)
			self.k4_status.acc = k4_force.fdiv(self.k4_status.mass)


			new_status_obj = Ship.copy(last_status_obj)
			# xf = x + (dt/6.0)*(v1 + 2*v2 + 2*v3 + v4)
			new_status_obj.pos = last_status_obj.pos.add(
				self.k1_status.vel.add(
					self.k2_status.vel.fmul(2.0)).add(
					self.k3_status.vel.fmul(2.0)).add(
					self.k4_status.vel).fmul(t_step/6.0))
			# vf = v + (dt/6.0)*(a1 + 2*a2 + 2*a3 + a4)
			new_status_obj.vel = last_status_obj.vel.add(
				self.k1_status.acc.add(
					self.k2_status.acc.fmul(2.0)).add(
					self.k3_status.acc.fmul(2.0)).add(
					self.k4_status.acc).fmul(t_step/6.0))

			#progress.update(t_curr)
			
			i += 1
			t_curr = t_step * i #instead of adding, to increase accuracy
			
			orbit_path.put_items(t_curr, t_step, new_status_obj)
			last_status_obj = new_status_obj


def vector2d_perpendicular(vec, i):
	"""i is either 0 or 1 depending on which perpendicular you want to select"""
	if i == 0:
		return Vector3(-vec.y, vec.x, vec.z)
	else:
		return Vector3(vec.y, -vec.x, vec.z)

class Star(object):
	def __init__(self, pos, color):
		self.pos = pos
		self.color = color

	def twinkle(self, dt):
		pass

class Viewport(Widget):
	def __init__(self, **kwargs):
		super(Viewport, self).__init__()
		self.render_objects = []
		if 'render_objects' in kwargs:
			self.render_objects += kwargs['render_objects']

		if 'world_scale' in kwargs:
			self.world_scale = kwargs['world_scale']
		else:
			self.world_scale = 1

		self.ship_scale = 15
		self.num_stars = 100
		self.star_brightness = 0.8
		self.first_update = True

	def add_render_object(self, obj):
		self.render_objects.append(obj)

	def remove_render_object(self, obj):
		self.render_objects.remove(obj)

	def update(self, dt):
		if self.first_update:
			self.stars = []
			i = 0
			while i < self.num_stars:
				pos = (random.randrange(self.parent.width), 
					random.randrange(self.parent.height))
				color = Color(0, 0, random.random()*self.star_brightness,mode='hsv')
				self.stars.append(Star(pos, color))
				i+=1
			self.first_update = False
		
		self.canvas.clear()
		self.pos = self.parent.pos
		with self.canvas:
			center_vec = Vector3(self.parent.center_x, self.parent.center_y, 0)

			for star in self.stars:
				Color(star.color.h, star.color.s, star.color.v, mode='hsv')
				Point(points=star.pos)

			for obj in self.render_objects:
				#Render Planets
				Color(1,1,1)
				if type(obj) == Planet:
					rad = obj.radius*self.world_scale
					Ellipse(pos=(self.parent.center_x-rad, self.parent.center_y-rad),
						size=(rad*2, rad*2))
					# Ellipse(pos=self.parent.center, size=(10,10))

				#Render Orbit Path
				if obj.__class__ == OrbitPath:
					Color(1,0,0)
					point_list = obj.get_point_list(center_vec, self.world_scale, step=2)
					Line(points=point_list, dash_length=2, dash_offset=2, width=1.2)

				#Render the Ship
				if type(obj) == Ship:
					#Ship Body
					Color(0, 1, 0)
					ship_pos = center_vec.add(obj.pos.fmul(self.world_scale))
					ship_dir = obj.direction
					ship_dir_left = vector2d_perpendicular(ship_dir, 0)
					ship_dir_right = vector2d_perpendicular(ship_dir, 1)
					ship_tri_front = ship_pos.add(ship_dir.fmul(self.ship_scale))
					ship_tri_back = ship_pos.add(ship_dir.fmul(-0.5* self.ship_scale))
					ship_tri_back_left = ship_tri_back.add(ship_dir_left.fmul(self.ship_scale*0.4))
					ship_tri_back_right = ship_tri_back.add(ship_dir_right.fmul(self.ship_scale*0.4))

					Triangle(points=(ship_tri_front.x, ship_tri_front.y,
						ship_tri_back_left.x, ship_tri_back_left.y,
						ship_tri_back_right.x, ship_tri_back_right.y))

					#Engine Fire
					Color(1, 0.4, 0.4)
					fire_size_x = self.ship_scale*0.3
					fire_size_y = self.ship_scale*1e-1*obj.get_effective_engine_power()
					fire_pos_top_left = ship_tri_back.add(ship_dir_left.fmul(fire_size_x/2.0))
					fire_pos_top_right = ship_tri_back.add(ship_dir_left.fmul(-fire_size_x/2.0))
					fire_pos_bottom = ship_tri_back.add(ship_dir.fmul(-fire_size_y))
					Triangle(points=(fire_pos_top_left.x, fire_pos_top_left.y,
						fire_pos_bottom.x, fire_pos_bottom.y,
						fire_pos_top_right.x, fire_pos_top_right.y))

					#Ship Engine
					Color(0.4, 0.4, 0.4)
					engine_size_x = self.ship_scale*0.5
					engine_size_y = self.ship_scale*0.5
					engine_pos_top_left = ship_tri_back.add(ship_dir_left.fmul(engine_size_x/2.0))
					engine_pos_top_right = ship_tri_back.add(ship_dir_left.fmul(-engine_size_x/2.0))
					engine_pos_bottom_left = engine_pos_top_left.add(ship_dir.fmul(-engine_size_y/2.0))
					engine_pos_bottom_right = engine_pos_top_right.add(ship_dir.fmul(-engine_size_y/2.0))

					Triangle(points=(engine_pos_top_left.x, engine_pos_top_left.y,
						engine_pos_top_right.x, engine_pos_top_right.y,
						engine_pos_bottom_right.x, engine_pos_bottom_right.y))
					Triangle(points=(engine_pos_top_left.x, engine_pos_top_left.y,
						engine_pos_bottom_left.x, engine_pos_bottom_left.y,
						engine_pos_bottom_right.x, engine_pos_bottom_right.y))


def to_label_str(number):
	return "{0:.2f}".format(number)			

class KivyOrbiter(Widget):
	def __init__(self):
		super(KivyOrbiter, self).__init__()

		self.ship_grab = False
		self.ship_grab_vec = Vector3()
		self.ship_grab_init_pos = Vector3()

		self.update_orbit_flag = True

		self.preview_thread = None
		self.orbit_t = 0.0
		self.physics_dt = 1/30.0
		
		self.planet = Planet(vel=Vector3(),
			pos=Vector3(),
			acc=Vector3(), 
			mass=5.972E25,
			radius=6371e3,
			show_vel=False,
			color='white')
		self.ship = Ship(direction=Vector3.Y(), 
						vel=Vector3(0, 7.66e3*3.5, 0),
                        pos=Vector3(float(6371e3 + 800e3), 0, 0),
                        acc=Vector3(),
                        mass=450000.0,
                        radius=0.0,
                        color="green",
                        show_vel=True)

		self.ship_orbit = RK4PreviewOrbit(self.planet, self.ship)
		self.ship_orbit_path = OrbitPath(self.ship)
		self.viewport = Viewport(render_objects=[self.planet, self.ship_orbit_path, self.ship], world_scale=5e-6)

		self.fuel_label = Label(text=str(to_label_str(self.ship.fuel)), color=[1,1,1,1])
		self.add_widget(self.viewport, 0)
		self.add_widget(self.fuel_label, 1)
		Clock.schedule_interval(self.update, 1.0/30.0)

	def on_touch_down(self, touch):
		if not self.ship_grab:
			touch.grab(self)
			self.ship_grab = True
			self.ship_grab_init_pos = Vector3(touch.spos[0], touch.spos[1], 0)

	def on_touch_move(self, touch):
		if touch.grab_current is self:
			grab_pos = Vector3(touch.spos[0], touch.spos[1], 0)
			self.ship_grab_vec = grab_pos.sub(self.ship_grab_init_pos)

	def on_touch_up(self, touch):
		if touch.grab_current is self:
			touch.ungrab(self)
			self.ship_grab = False
			self.ship_grab_vec = Vector3()
			self.update_orbit_flag = True
				

	def update_orbit(self, physics_dt):
		orbit_max_t = physics_dt*600000
		self.orbit_t = 0.0
		self.ship_orbit.triggger_calc_orbit(orbit_max_t, 
			physics_dt*3000, self.ship_orbit_path, reset=True)
		
		
	def update(self, dt):
		
		if self.update_orbit_flag:
			self.update_orbit(self.physics_dt)
			self.update_orbit_flag = False

		#wait until orbit path has been populated
		while len(self.ship_orbit_path) < 2:
			pass

		status_obj = self.ship_orbit_path.get_status_interp(self.orbit_t, trim=False)
		if status_obj is not None:
			self.ship.pos = status_obj.pos
			self.ship.vel = status_obj.vel

		grab_mag = self.ship_grab_vec.mag()
		if grab_mag > 0.001:
			self.ship.direction = self.ship_grab_vec.norm()
			self.ship.engine_power = grab_mag * 100
			self.update_orbit(self.physics_dt)
		else:
			self.ship.engine_power = 0.0

		self.ship.update(dt)
		self.fuel_label.text = str(to_label_str(self.ship.fuel))
		
		self.viewport.update(dt)

		#if self.orbit_t > 1000:
		#	self.app.stop()

		self.orbit_t += dt*400



class KivyOrbiterApp(App):
	def build(self):
		return KivyOrbiter()


if __name__ == '__main__':
	Config.set('graphics', 'width', '800')
	Config.set('graphics', 'height', '480')

	#import yappi
	#yappi.start()
	KivyOrbiterApp().run()
	#yappi.get_func_stats().print_all()
	#yappi.get_thread_stats().print_all()
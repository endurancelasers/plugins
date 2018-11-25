#Version controll: last edited by 14.04.2016 17:44

#!/usr/bin/env python 

from __future__ import division
from math import sqrt, cos, sin, tan, atan, atan2, pi, hypot, copysign, pow, ceil
import inkex, simplestyle, simplepath
import cubicsuperpath, simpletransform, bezmisc

from subprocess import Popen, PIPE

import os
import sys
#import gettext
#_ = gettext.gettext

 
### Check if inkex has errormsg (0.46 version doesnot have one.) Could be removed later.
#if "errormsg" not in dir(inkex):
#    inkex.errormsg = lambda msg: sys.stderr.write((unicode(msg) + "\n").encode("UTF-8"))

inf = float("inf")
pi2 = pi*2
def print_(s):
    f = open("/home/dingo/gcodes/hatching/hatch.log", "a+")
    f.write("%s\n" % s)
    f.close()

def lerp(a, b, t):
    return a+t*(b-a)
    
def ilerp(a, b, v):
    return (v-a)/(b-a)
    
def lerp_a(p0, p1, t):
    return [lerp(p0[0], p1[0], t), lerp(p0[1], p1[1], t)]
        
def dot(p0, p1):
    return p0[0]*p1[0]+p0[1]*p1[1]
    
def Arg(dx, dy):
    a = 0
    if dx == 0:
	    a = pi/2
    elif dx>0:
        a = atan(abs(dy/dx))
    else:
        a = pi-atan(abs(dy/dx))

    if dy>=0:
        return a
    else:
        return -a
        
def arc_center(p1, p2, p3):
    x1 = p1[0]
    y1 = p1[1]

    x2 = p2[0]
    y2 = p2[1]
    
    x3 = p3[0]
    y3 = p3[1]

    x12 = x1-x2
    y12 = y1-y2
    
    x23 = x2-x3
    y23 = y2-y3
    
    x31 = x3-x1
    y31 = y3-y1
    
    z1 = x1*x1+y1*y1
    z2 = x2*x2+y2*y2
    z3 = x3*x3+y3*y3
    
    zx = y12*z3+y23*z1+y31*z2
    zy = x12*z3+x23*z1+x31*z2
    z = x12*y31-y12*x31
    
    try:
        cx = -zx/2/z
        cy = zy/2/z
        
        x3c = x3 - cx
        y3c = y3 - cy
        
        p = x3c*y3c - y12*x3c
        
        return [cx, cy, copysign(hypot(x1-cx, y1-cy), -z)]
#        return [cx, cy, copysign(hypot(x1-cx, y1-cy), -z if cx*cy>=0 else z)]
    except ZeroDivisionError:
        return []
        
def px_to_mm(val):
    return val / 3.5433070660

def mm_to_px(val):
    return val * 3.5433070660
    
def frange(start, stop, step):
    cnt=start
    while cnt<=stop:
        yield cnt
        cnt+=step

class Polynome:
    def __init__(self, coeffs):
        self.coeffs = coeffs[:]
        
#        print_(str(self.coeffs))
        ok = False
        while len(self.coeffs) > 0 and not ok:
            if abs(self.coeffs[-1]) < 0.01:
                self.coeffs = self.coeffs[:-1]
                ok = False
            else:
                ok = True
        
    def value(self, t):
        _t = t
        value = self.coeffs[0]
        
        if len(self.coeffs)>1:
            for c in self.coeffs[1:]:
                value += _t*c
                _t *= t
        
        return value
        
    def derivative(self, t):
        if len(self.coeffs) < 2:
            return 0
            
        _t = t
        derivative = self.coeffs[1]
        factor = 2
        
        if len(self.coeffs)>2:
            for c in self.coeffs[2:]:
                derivative += factor*_t*c
                _t *= t
                factor+=1
        
        return derivative
        
    def derivative2(self, t):
        if len(self.coeffs) < 3:
            return 0
            
        _t = t
        factor = 2
        derivative2 = factor*self.coeffs[2]
        
        if len(self.coeffs) > 3:
            for c in self.coeffs[3:]:
                derivative2 += factor*(factor+1)*_t*c
                _t *= t
                factor+=1
        
        return derivative2
        
    def solve(self):
        result = []
        
        def solveLinear(a, b):
            result = []
#            print_("Solving linear")
            if a !=0:
                result = [-b/a]
            
            return result
        
        def solveQuadratic(a, b, c):
            result = []
#            print_("Solving quadratic")
            b /= a
            c /= a
            
            d = b**2-4*c
            
            if d>0:
                d = sqrt(d)
                result = [(-b-d)/2, (-b+d)/2]
            elif d == 0:
                result = [-b/2]
            
            return result
            
        def solveCubic(a, b, c, d):
            result = []
#            print_("Solving cubic")
            b /= a
            c /= a
            d /= a
            
            p = -b*b/3+c
            q = 2*b*b*b/27-b*c/3+d
            
            d1 = q*q/4+p*p*p/27
            
            if d1 >= 0:
                d1 = sqrt(d1)
                u = -q/2+d1
                u = copysign(abs(u)**(1/3), u)
                v = -q/2-d1
                v = copysign(abs(v)**(1/3), v)
                X = u+v-b/3
            else:
                X = 2*sqrt(-p/3)*cos(Arg(-q/2, sqrt(abs(d1)))/3)-b/3
                
            result = solveQuadratic(1, X+b, X*X+b*X+c)
            
            f = True
            for r in result:
                if abs(r - X)<0.001:
                    f = False
                    
            if f:
                result.append(X)
            
            return result

        if len(self.coeffs) > 3 and self.coeffs[3] != 0:
            result = solveCubic(self.coeffs[3], self.coeffs[2], self.coeffs[1], self.coeffs[0])
        elif len(self.coeffs) > 2 and self.coeffs[2] != 0:
            result = solveQuadratic(self.coeffs[2], self.coeffs[1], self.coeffs[0])
        elif len(self.coeffs) > 1 and self.coeffs[1] != 0:
            result = solveLinear(self.coeffs[1], self.coeffs[0])
        elif len(self.coeffs) > 0:
            result = [self.coeffs[0]]

        return result
                
class Bezier2d:
    def __init__(self, pts, linearity_tolerance = 1, length_threshold = 0.001, simplify = True):
        self.pts = pts[:]
        
        if len(self.pts) >2:
            if self.pts[0][0] == self.pts[1][0] and self.pts[0][1]==self.pts[1][1]:
                self.pts = self.pts[1:]
                
            if self.pts[-2][0] == self.pts[-1][0] and self.pts[-2][1]==self.pts[-1][1]:
                self.pts = self.pts[:-1]
    
        self.p_x = Polynome(self.calc_coeffs([p[0] for p in self.pts]))
        self.p_y = Polynome(self.calc_coeffs([p[1] for p in self.pts]))
    
#        if len(self.pts) > 2 and abs(self.curvature(0.5)+self.curvature(0)+self.curvature(1))*hypot(self.pts[0][0]-self.pts[-1][0], self.pts[0][1]-self.pts[-1][1])/3 < linearity_tolerance:
#        if len(self.pts) > 2 and abs(self.pts[0][0]*self.pts[1][1]-self.pts[0][1]*self.pts[1][0]) < linearity_tolerance and abs(self.pts[-1][0]*self.pts[-2][1]-self.pts[-1][1]*self.pts[-2][0]) < linearity_tolerance:
#        if simplify and len(self.pts)>2:
            
#            ends_dir = atan2(self.pts[-1][0]-self.pts[0][0], self.pts[-1][1]-self.pts[0][1])
#            sl_dir = atan2(self.pts[1][0]-self.pts[0][0], self.pts[1][1]-self.pts[0][1])
#            el_dir = atan2(self.pts[-1][0]-self.pts[-2][0], self.pts[-1][1]-self.pts[-2][1])
            
#            s_angle = abs(ends_dir - sl_dir)*180/pi
#            e_angle = abs(el_dir - ends_dir)*180/pi
           
#            if (s_angle < linearity_tolerance or s_angle > 180-linearity_tolerance) and (e_angle < linearity_tolerance or e_angle > 180-linearity_tolerance):
#                self.pts = [self.pts[0], self.pts[-1]]
            
        if simplify and len(self.pts) == 2:
            if px_to_mm(hypot(self.pts[0][0]-self.pts[-1][0], self.pts[0][1]-self.pts[-1][1])) >= length_threshold:
                self.p_x = Polynome(self.calc_coeffs([p[0] for p in self.pts]))
                self.p_y = Polynome(self.calc_coeffs([p[1] for p in self.pts]))
            else:
                self.p_x = None
                self.p_y = None
                self.pts = []
                
        if len(self.pts)>0:
            self.bbox = {"left":inf, "right":-inf, "top":-inf, "bottom":inf}
            
            for pt in self.pts:
                self.bbox["left"] = min(self.bbox["left"], pt[0])
                self.bbox["right"] = max(self.bbox["right"], pt[0])
                self.bbox["top"] = max(self.bbox["top"], pt[1])
                self.bbox["bottom"] = min(self.bbox["bottom"], pt[1])
        else:
            self.bbox = {"left":0, "right":0, "top":0, "bottom":0}
            
    def x(self, t):
        return self.p_x.value(t)

    def y(self, t):
        return self.p_y.value(t)
        
    def point(self, t):
        return [self.x(t), self.y(t)]
        
    def dx(self, t):
        return self.p_x.derivative(t)

    def d2x(self, t):
        return self.p_x.derivative2(t)

    def dy(self, t):
        return self.p_y.derivative(t)

    def d2y(self, t):
        return self.p_y.derivative2(t)
    
    def curvature(self, t):
        dx = self.dx(t)
        dy = self.dy(t)
        
        cl = hypot(dx, dy)**3
        
        t2 = t*t
        
        d2x = self.d2x(t)
        d2y = self.d2y(t)
        
        try:
            return (d2y*dx - d2x*dy)/cl
        except ZeroDivisionError:
            return 1e100
        
    def curv_radius(self, t):
        dx = self.dx(t)
        dy = self.dy(t)
        
        d2x = self.d2x(t)
        d2y = self.d2y(t)
        
        vp = d2y*dx - d2x*dy
        
        cl = hypot(dx, dy)**3
        
        try:
            return cl/vp
        except ZeroDivisionError:
            return 1e100
        
    def slope(self, t):
        return atan2(self.dy(t), self.dx(t))
    
    def curv_center_rel(self, t):
        dx = self.dx(t)
        dy = self.dy(t)
        
        d2x = self.d2x(t)
        d2y = self.d2y(t)
        
        t2 = t*t
        
        nom = dx*dx+dy*dy
        # denom = dx*d2y*d2y - d2x*d2x*dy
        denom = dx*self.d2y(t2) - self.d2x(t2)*dy
        
        try:
            return [-dy*nom/denom, dx*nom/denom]
        except ZeroDivisionError:
            return []

    def curv_center(self, t):
        rel_center = self.curv_center_rel(t)
        if len(rel_center)>0:
            return [self.x(t)+rel_center[0], self.y(t)+rel_center[1]]
        else:
            return []
            
    def dist_to_arc(self, cx, cy, radius):
        max_dist = 0
        for i in xrange(11):
            t = float(i)/10
            x, y = self.point(t)
            dist = abs(hypot(x-cx, y-cy) - radius)
            
            max_dist = max(max_dist, dist)
        
        return max_dist
        
    def dist_to_line(self, p1, p2):
        denom = hypot(p2[0] - p1[0], p2[1] - p1[1])
        
        if denom == 0:
            return 0
        
        y01 = p1[1] - p2[1]
        x10 = p2[0] - p1[0]
        a = p1[0]*p2[1] - p2[0]*p1[1]
        
        max_dist = 0
        for i in xrange(11):
            t = float(i)/10
            x, y = self.point(t)
            dist = (x*y01+y*x10+a)/denom
            
            max_dist = max(max_dist, dist)
            
        return max_dist
        
    def intersect_to_line(self, a, b):
        if a[0] == b[0] and a[1] == b[1]:
            return []
    
        result = []
        
        def intersect_line_to_line(p0, p1, a, b):
            result = []
            
            ua_t=(b[0]-a[0])*(p0[1]-a[1])-(b[1]-a[1])*(p0[0]-a[0])
            ub_t=(p1[0]-p0[0])*(p0[1]-a[1])-(p1[1]-p0[1])*(p0[0]-a[0])
            u_b=(b[1]-a[1])*(p1[0]-p0[0])-(b[0]-a[0])*(p1[1]-p0[1])
            
            if u_b != 0:
                ua=ua_t/u_b
                ub=ub_t/u_b
                
                if 0<=ua and ua<=1 and 0<=ub and ub<=1:
                    result = [ub]
            
            return result
            
        def intersect_quad_to_line(p0, p1, p2, a, b):
            result = []
            
            _min = [min(a[0], b[0]), min(a[1], b[1])]
            _max = [max(a[0], b[0]), max(a[1], b[1])]
            
            _a = [-2*p1[0], -2*p1[1]]
            c2 = [p0[0]+_a[0]+p2[0], p0[1]+_a[1]+p2[1]]
            
            _a = [-2*p0[0], -2*p0[1]]
            _b = [2*p1[0], 2*p1[1]]
            c1 = [_a[0]+_b[0], _a[1]+_b[1]]
            
            c0 = p0[:]
            n = [a[1]-b[1], b[0]-a[0]]
            
            cl = a[0]*b[1]-b[0]*a[1]
            poly = Polynome([dot(c0, n)+cl, dot(c1, n), dot(c2, n)])
            
            roots = poly.solve()
            
            for t in roots:
                if 0<=t and t<=1:
                    lt = 0
                    
                    if abs(a[0]-b[0]) > abs(a[1]-b[1]):
                        lt = ilerp(a[0], b[0], self.x(t))
                    else:
                        lt = ilerp(a[1], b[1], self.y(t))
                        
#                    if 0<=lt and lt <=1:    
                    if True:
                        result.append(lt)
                        
            return result
            
        def intersect_cubic_to_line(p0, p1, p2, p3, a, b):    
            result = []
            
            _min = [min(a[0], b[0]), min(a[1], b[1])]
            _max = [max(a[0], b[0]), max(a[1], b[1])]
            
            _a = [-1*p0[0], -1*p0[1]]
            _b = [3*p1[0], 3*p1[1]]
            _c = [-3*p2[0], -3*p2[1]]
            _d = [_a[0]+_b[0]+_c[0]+p3[0], _a[1]+_b[1]+_c[1]+p3[1]]
            
            c3 = _d[:]
            
            _a = [3*p0[0], 3*p0[1]]
            _b = [-6*p1[0], -6*p1[1]]
            _c = [3*p2[0], 3*p2[1]]
            _d = [_a[0]+_b[0]+_c[0], _a[1]+_b[1]+_c[1]]
            
            c2 = _d[:]
            
            _a = [-3*p0[0], -3*p0[1]]
            _b = [3*p1[0], 3*p1[1]]
            _c = [_a[0]+_b[0], _a[1]+_b[1]]
            
            c1 = _c[:]
            c0 = p0[:]
            
            n = [a[1]-b[1], b[0]-a[0]]
            cl = a[0]*b[1]-b[0]*a[1]
            
            poly = Polynome([dot(c0, n)+cl, dot(c1, n), dot(c2, n), dot(c3, n)])
            
#            sys.stderr.write("%s\n" % str(poly.coeffs))
            
            roots = poly.solve()

            for t in roots:
                if 0<=t and t<=1:
                    lt = 0
                    
                    if abs(a[0]-b[0]) > abs(a[1]-b[1]):
                        lt = ilerp(a[0], b[0], self.x(t))
                    else:
                        lt = ilerp(a[1], b[1], self.y(t))
                    
#                    if 0<=lt and lt<=1:    
                    if True:
                        result.append(lt)

            return result
            
        def line_intersects_to_bbox(left, right, top, bottom, a, b):
            result = False
            
            if len(intersect_line_to_line([left, top], [right, top], a, b)) > 0:
                result = True
            elif len(intersect_line_to_line([right, top], [right, bottom], a, b)) > 0:
                result = True
            elif len(intersect_line_to_line([left, bottom], [right, bottom], a, b)) > 0:
                result = True
            elif len(intersect_line_to_line([left, top], [left, bottom], a, b)) > 0:
                result = True
            
            return result
    
#        if line_intersects_to_bbox(self.bbox["left"], self.bbox["right"], self.bbox["top"], self.bbox["bottom"], a, b):
        if True:
#            sys.stderr.write("Intersects to BBox\n")
            if len(self.pts)>3:
                result = intersect_cubic_to_line(self.pts[0], self.pts[1], self.pts[2], self.pts[3], a, b)
            elif len(self.pts)>2:
                result = intersect_quad_to_line(self.pts[0], self.pts[1], self.pts[2], a, b)
            elif len(self.pts)>1:
                result = intersect_line_to_line(self.pts[0], self.pts[1], a, b)
            
        return result
        
    def approximate(self, tolerance = 0.0001, linearity_tolerance = 1):
        buf = []
        
        if len(self.pts) == 0:
            return []
            
        if len(self.pts) == 2:
            return [["line", self.pts[0], self.pts[1]]]
        else:
            if self.pts[0][0] == self.pts[-1][0] and self.pts[0][1] == self.pts[-1][1]:
                m1 = self.point(0.2)
                m2 = self.point(0.8)
            else:
                m1 = self.pts[0][:]
                m2 = self.pts[-1][:]
            
            m3 = self.point(0.4)
            m4 = self.point(0.6)
            center = arc_center(m1, [(m3[0]+m4[0])/2, (m3[1]+m4[1])/2], m2)
                        
            if len(center) > 0:
                cx, cy, r = center

                start_angle = atan2(m1[1]-cy, m1[0]-cx)
                end_angle = atan2(m2[1]- cy, m2[0] - cx)
                angle = (end_angle - start_angle) % pi2
                
                if abs(r) < 10000 and angle*180/pi > max(0.1,linearity_tolerance):
                
                    r1 = hypot(self.pts[0][0]-cx, self.pts[0][1]-cy)
                    r2 = hypot(self.pts[-1][0]-cx, self.pts[-1][1]-cy)                    
            
                    if abs(r1-r2) > 0.0001 or self.dist_to_arc(cx, cy, r1) > tolerance:
                        f, s = self.split(0.5)
                        buf+=f.approximate(tolerance, linearity_tolerance)
                        buf+=s.approximate(tolerance, linearity_tolerance)
                        
                        return buf
                    else:
                        return [["arc", cx-self.pts[0][0], cy-self.pts[0][1], copysign(r1,r), self.pts[0], self.pts[-1]]]
                else:
                    return [["line", self.pts[0], self.pts[-1]]]
            else:
                return [["line", self.pts[0], self.pts[-1]]]
                
    def split(self, t):
        if len(self.pts) == 0:
            return []
    
        split_x = self.x(t)
        split_y = self.y(t)
            
        if len(self.pts) == 2:
            return [Bezier2d([self.pts[0], [split_x, split_y]]), Bezier2d([[split_x, split_y], self.pts[1]])]
        elif len(self.pts) == 3:
            x01 = lerp(self.pts[0][0], self.pts[1][0], t)
            y01 = lerp(self.pts[0][1], self.pts[1][1], t)
            
            x12 = lerp(self.pts[1][0], self.pts[2][0], t)
            y12 = lerp(self.pts[1][1], self.pts[2][1], t)
            
            return [Bezier2d([self.pts[0], [x01, y01], [split_x, split_y]]), Bezier2d([[split_x, split_y], [x12, y12], self.pts[2]])]
        elif len(self.pts) == 4:
            x01 = lerp(self.pts[0][0], self.pts[1][0], t)
            y01 = lerp(self.pts[0][1], self.pts[1][1], t)
            
            x12 = lerp(self.pts[1][0], self.pts[2][0], t)
            y12 = lerp(self.pts[1][1], self.pts[2][1], t)
            
            x23 = lerp(self.pts[2][0], self.pts[3][0], t)
            y23 = lerp(self.pts[2][1], self.pts[3][1], t)
            
            x012 = lerp(x01, x12, t)
            y012 = lerp(y01, y12, t)
            
            x123 = lerp(x12, x23, t)
            y123 = lerp(y12, y23, t)
            
            return [Bezier2d([self.pts[0], [x01, y01], [x012, y012], [split_x, split_y]]), Bezier2d([[split_x, split_y], [x123, y123], [x23, y23], self.pts[3]])]
    
    def calc_line_coeffs(self, p0, p1):
        return [p0, p1-p0]
    
    def calc_quad_coeffs(self, p0, p1, p2):
        return [p0, 2*(p1-p0), p0+p2-2*p1]
    
    def calc_cubic_coeffs(self, p0, p1, p2, p3):
        return [p0, 3*(p1-p0), 3*(p0+p2-2*p1), 3*(p1-p2)+p3-p0]
    
    def calc_coeffs(self, pts):
        if len(pts) == 0:
            return
            
        if len(pts) == 2:
            return self.calc_line_coeffs(*pts)
        elif len(pts) == 3:
            return self.calc_quad_coeffs(*pts)
        elif len(pts) == 4:
            return self.calc_cubic_coeffs(*pts)

import re
class PathParser:
    def __init__(self, path_string):
        COMMANDS = set('MmZzLlHhVvCcSsQqTtAa')
        COMMAND_RE = re.compile("([MmZzLlHhVvCcSsQqTtAa])")
        FLOAT_RE = re.compile("[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?")
        
        self.cmds = []
        
        for x in COMMAND_RE.split(path_string):
            cmd = []
            if x in COMMANDS:
                self.cmds.append(x)

            for token in FLOAT_RE.findall(x):
                cmd.append(token)
                    
            if len(cmd) > 0:
                self.cmds.append(cmd[:])
            
    def print_cmds(self):
        sys.stderr.write(str(self.cmds))

class GcodeGenerator():
    def __init__(self, laser_on="M106", laser_off="M107", laser_power = 254, laser_on_pause=100, laser_off_pause = 100, feed=100, preambule = "", postambule = "", path_preambule = "", path_postambule = "", delete_comments = ""):
        self.laser_on = laser_on
        self.laser_off = laser_off
        self.laser_power = laser_power
        self.laser_on_pause = laser_on_pause
        self.laser_off_pause = laser_off_pause
        self.feed = feed
        self.preambule = preambule
        self.postambule = postambule
        self.path_preambule = path_preambule
        self.path_postambule = path_postambule
	self.delete_comments = delete_comments
	
        self.buf = ""
        
    def add_path(self, appr, name=""):
        if not appr:
            return
        
        start_detected = False
        path_start = []
        
        buf = ""
        for path in appr:
            if path[0] == "line":
                typ, start, end = path
                buf += "G01 X%.4f Y%.4f F%d\n" % (px_to_mm(end[0]), px_to_mm(self.height-end[1]), self.feed)
            elif path[0] == "arc":
                typ, cx, cy, r, start, end = path
                buf += "G%s X%.4f Y%.4f I%.4f J%.4f F%d\n" % ("03" if r<0 else "02", px_to_mm(end[0]), px_to_mm(self.height-end[1]), px_to_mm(cx), px_to_mm(-cy), self.feed)

            if not start_detected:
                start_detected = True
                path_start = start[:]
               
       
	if self.delete_comments:
		start_seq = "G00 X%.4f Y%.4f\nG04 P%d\n%s S%d\n%s\n%s\n" % (px_to_mm(path_start[0]), px_to_mm(self.height-path_start[1]), self.laser_on_pause, self.laser_on, self.laser_power, self.path_preambule, " "+name if name else "")
		end_seq = "G04 P%d\n%s S0\n%s\n" % (self.laser_off_pause, self.laser_off, self.path_postambule)
	else: 
		start_seq = "G00 X%.4f Y%.4f\nG04 P%d\n%s S%d\n%s\n;(Path%s start)\n" % (px_to_mm(path_start[0]), px_to_mm(self.height-path_start[1]), self.laser_on_pause, self.laser_on, self.laser_power, self.path_preambule, " "+name if name else "")
		end_seq = "G04 P%d\n%s S0\n%s\n;(Path end)\n\n" % (self.laser_off_pause, self.laser_off, self.path_postambule)
            
        self.buf += start_seq + buf + end_seq
        
    def start(self, height):
        self.height = float(height[:-2])
        if height[-2:] == "mm":
            self.height = mm_to_px(self.height)
        
    def get_preambule(self):
		if self.delete_comments:
			return "\n%s\nG00\nG04 P%d\n%s S0\nG90 G21\n" % (self.preambule, self.laser_off_pause, self.laser_off)
		else: 
			return ";Endurance Robots Laser Lab. http://endurancerobots.com\n;(G-Code file)\n%s\nG00\nG04 P%d\n%s S0\nG90 G21\n" % (self.preambule, self.laser_off_pause, self.laser_off)
			
    def get_gcode(self):
        return self.buf

    def get_postambule(self):
		if self.delete_comments:
			return "G04 P%d\n%s S0\nG00\n%s\nM02\n\n" % (self.laser_off_pause, self.laser_off, self.postambule)
		else: 
			return "G04 P%d\n%s S0\nG00\n%s\nM02\n\n;Endurance Robots Laser Lab. http://endurancerobots.com" % (self.laser_off_pause, self.laser_off, self.postambule)
			
    def end(self):
        return "%s\n%s\n%s" % (self.get_preambule(), self.get_gcode(), self.get_postambule())

    def gen_svg(self, path):
        start_detected = False
        seg_start = []
        d_buf = ""
        
        for seg in path:
            if seg[0] == "line":
                typ, start, end = seg
                d_buf += "L %f,%f " % (end[0], end[1])
            elif seg[0] == "arc":
                typ, cx, cy, r, start, end = seg
                radius = hypot(cx, cy)
                d_buf += "A %f,%f 0 0 %d %f,%f " %(radius, radius, 1 if (r > 0 and r < 1e100) else 0, end[0], end[1])
            if not start_detected:
                start_detected = True
                seg_start = start
        
        if len(seg_start) > 0:        
            d_buf = "M %f,%f %s " % (seg_start[0], seg_start[1], d_buf)
        
        return d_buf

class RGcodeEffect(inkex.Effect):
    def __init__(self):
        inkex.Effect.__init__(self)
		
	self.OptionParser.add_option("--delete-comments", 
						action="store", 
						type="inkbool",    
						dest="delete_comments", 
						default=True,
						help="Delete all comments")
		
        self.OptionParser.add_option("","--add-numeric-suffix-to-filename", 
						action="store", 
						type="inkbool",    
						dest="add_numeric_suffix_to_filename", 
						default=True,
						help="Add numeric suffix to filename")
		
        self.OptionParser.add_option("--directory",
                        action="store", type="string", 
                        dest="directory", default="",
                        help="directory")

        self.OptionParser.add_option("--file",
                        action="store", type="string", 
                        dest="file", default="",
                        help="file")
        
        self.OptionParser.add_option("-o", "--laser_on",
                        action="store", type="string", 
                        dest="laser_on", default="M106",
                        help="laser ON command")

        self.OptionParser.add_option("-x", "--laser_off",
                        action="store", type="string", 
                        dest="laser_off", default="M107",
                        help="laser OFF command")
                        
        self.OptionParser.add_option("-s", "--job_preambule",
                        action="store", type="string", 
                        dest="job_preambule", default="",
                        help="job preambule")

        self.OptionParser.add_option("-d", "--job_postambule",
                        action="store", type="string", 
                        dest="job_postambule", default="",
                        help="job postambule")

        self.OptionParser.add_option("-b", "--path_preambule",
                        action="store", type="string", 
                        dest="path_preambule", default="",
                        help="path preambule")

        self.OptionParser.add_option("-e", "--path_postambule",
                        action="store", type="string", 
                        dest="path_postambule", default="",
                        help="path postambule")
                        
        self.OptionParser.add_option("-p", "--laser_power",
                        action="store", type="int", 
                        dest="laser_power", default=254,
                        help="laser power")
        
        self.OptionParser.add_option("", "--laser_on_pause",
                        action="store", type="int", 
                        dest="laser_on_pause", default=100,
                        help="pause before path start")

        self.OptionParser.add_option("", "--laser_off_pause",
                        action="store", type="int", 
                        dest="laser_off_pause", default=100,
                        help="pause after path end")

        self.OptionParser.add_option("-f", "--feed_speed",
                        action="store", type="int", 
                        dest="feed_speed", default=100,
                        help="feed speed")
                        
        self.OptionParser.add_option("-a", "--approximation_tolerance",
                        action="store", type="float", 
                        dest="approximation_tolerance", default=400,
                        help="approximation tolerance")

        self.OptionParser.add_option("-l", "--linearity_tolerance",
                        action="store", type="float", 
                        dest="linearity_tolerance", default=5,
                        help="linearity tolerance")
                        
        self.OptionParser.add_option("-t", "--length_threshold",
                        action="store", type="float", 
                        dest="length_threshold", default=10,
                        help="length threshold")

        self.OptionParser.add_option("-m", "--show_markers",
                        action="store", type="inkbool",
                        dest="show_markers", default=True,
                        help="show markers")
                        
        self.OptionParser.add_option('', '--angle', type='float', dest='angle', default='45', help='Angle of lines', action="store")
        self.OptionParser.add_option('', '--distance', type='int', dest='distance', default='5', help='Distance between lines', action="store")
        self.OptionParser.add_option('', '--fill-method', type='string', dest='fill_method', default='lines', help='Fill method', action="store")
        self.OptionParser.add_option('', '--fill-rule', type='string', dest='rule', default='odd', help='Fill rule', action="store")
        
        self.OptionParser.add_option('', '--paths', type='inkbool', dest='paths', default=False, action="store")
	self.OptionParser.add_option('', '--filling', type='inkbool', dest='filling', default=False, action="store")
        self.OptionParser.add_option('', '--active-tab', type='string', dest='active_tab', default='', action="store")
       
		
    def export_gcode(self,gcode):
        gcode_pass = gcode

        f = open(os.path.join(self.options.directory,self.options.file), "w")
        f.write(gcode)
        f.close()
        
    def get_transforms(self,g):
        root = self.document.getroot()
        trans = []
        while (g!=root):
            if 'transform' in g.keys():
                t = g.get('transform')
                t = simpletransform.parseTransform(t)
                trans = simpletransform.composeTransform(t,trans) if trans != [] else t
#                print_(trans)
            g=g.getparent()
        return trans
        
    def apply_transforms(self,g,csp):
        trans = self.get_transforms(g)
        if trans != []:
            simpletransform.applyTransformToPath(trans, csp)
        return csp
        
    def get_info(self):
        self.selected_paths = {}
        self.paths = {}        
        self.orientation_points = {}
        self.layers = [self.document.getroot()]
        self.Zcoordinates = {}
        self.transform_matrix = {}
        self.transform_matrix_reverse = {}
        self.Zauto_scale = {}
        
        def recursive_search(g, layer, selected=False):
            items = g.getchildren()
            items.reverse()
            for i in items:
                if selected:
                    self.selected[i.get("id")] = i
                if i.tag == inkex.addNS("g",'svg') and i.get(inkex.addNS('groupmode','inkscape')) == 'layer':
                    self.layers += [i]
                    recursive_search(i,i)
#                elif i.get('gcodetools') == "Gcodetools orientation group" :
#                    points = self.get_orientation_points(i)
#                    if points != None :
#                        self.orientation_points[layer] = self.orientation_points[layer]+[points[:]] if layer in self.orientation_points else [points[:]]
#                        print_("Found orientation points in '%s' layer: %s" % (layer.get(inkex.addNS('label','inkscape')), points))
#                    else :
#                        self.error(_("Warning! Found bad orientation points in '%s' layer. Resulting Gcode could be corrupt!") % layer.get(inkex.addNS('label','inkscape')), "bad_orientation_points_in_some_layers")
                elif i.tag == inkex.addNS('path','svg'):
                    if "gcodetools"  not in i.keys() :
                        self.paths[layer] = self.paths[layer] + [i] if layer in self.paths else [i]  
                        if i.get("id") in self.selected :
                            self.selected_paths[layer] = self.selected_paths[layer] + [i] if layer in self.selected_paths else [i]  
                elif i.tag == inkex.addNS("g",'svg'):
                    recursive_search(i,layer, (i.get("id") in self.selected) )
                elif i.get("id") in self.selected :
                    pass
#                    self.error(_("This extension works with Paths and Dynamic Offsets and groups of them only! All other objects will be ignored!\nSolution 1: press Path->Object to path or Shift+Ctrl+C.\nSolution 2: Path->Dynamic offset or Ctrl+J.\nSolution 3: export all contours to PostScript level 2 (File->Save As->.ps) and File->Import this file."),"selection_contains_objects_that_are_not_paths")
                
                    
        recursive_search(self.document.getroot(),self.document.getroot())

    def get_defs(self) :
        self.defs = {}
        def recursive(g) :
            for i in g:
                if i.tag == inkex.addNS("defs","svg") : 
                    for j in i: 
                        self.defs[j.get("id")] = i
                if i.tag ==inkex.addNS("g",'svg') :
                    recursive(i)
        recursive(self.document.getroot())
    
    def get_bbox(self, ids):
        minX = float("inf")
        maxX = -float("inf")
        minY = float("inf")
        maxY = -float("inf")
        
        for id in ids:
            # query bounding box, UPPER LEFT corner (?)
    		q = {'x':0, 'y':0, 'width':0, 'height':0}
    		for query in q.keys():
    			p = Popen(
    				'inkscape --query-%s --query-id=%s "%s"' % (query, id, self.args[-1], ),
    				shell=True,
    				stdout=PIPE,
    				stderr=PIPE,
    				)
    			p.wait()
    			q[query] = p.stdout.read()
    
    		# get width, height, center of bounding box 
    		obj_width = float(q['width'])
    		obj_height = float(q['height'])
    		obj_x = float(q['x'])
    		obj_y = float(q['y'])
    		
    		minX = min(minX, min(obj_x, obj_x+obj_width))
    		maxX = max(maxX, max(obj_x, obj_x+obj_width))
    		minY = min(minY, min(obj_y, obj_y+obj_height))
    		maxY = max(maxY, max(obj_y, obj_y+obj_height))
    		
        return {"left": minX, "right": maxX, "top": maxY, "bottom": minY}
        
    def draw_path(self, path, stroke_width = 2, group = None):
        if self.options.show_markers:
            line_style = "fill:none;stroke:#FF0000; width:%dpx;marker-start:url(#DrawCurveMarker)" % stroke_width
            arc_style = "fill:none;stroke:#FF0000; width:%dpx;marker-start:url(#DrawCurveMarker)" % stroke_width
            arc_style_r = "fill:none;stroke:#FF0000; width:%dpx;marker-end:url(#DrawCurveMarker_r)" % stroke_width
        else:
            line_style = "fill:none;stroke:#FF0000; width:%dpx;marker-start:none" % stroke_width
            arc_style = "fill:none;stroke:#FF0000; width:%dpx;marker-start:none" % stroke_width
            arc_style_r = "fill:none;stroke:#FF0000; width:%dpx;marker-end:none" % stroke_width
        pi2 = pi*2
        
        self.get_defs()
        # Add marker to defs if it doesnot exists
        if "DrawCurveMarker" not in self.defs : 
            defs = inkex.etree.SubElement( self.document.getroot(), inkex.addNS("defs","svg"))
            marker = inkex.etree.SubElement( defs, inkex.addNS("marker","svg"), {"id":"DrawCurveMarker","orient":"auto","refX":"-8","refY":"-2.41063","style":"overflow:visible"})
            inkex.etree.SubElement( marker, inkex.addNS("path","svg"), 
                    {    "d":"m -6.55552,-2.41063 0,0 L -13.11104,0 c 1.0473,-1.42323 1.04126,-3.37047 0,-4.82126",
                        "style": "fill:#000044; fill-rule:evenodd;stroke-width:0.62500000;stroke-linejoin:round;"    }
                )
        if "DrawCurveMarker_r" not in self.defs : 
            defs = inkex.etree.SubElement( self.document.getroot(), inkex.addNS("defs","svg"))
            marker = inkex.etree.SubElement( defs, inkex.addNS("marker","svg"), {"id":"DrawCurveMarker_r","orient":"auto","refX":"8","refY":"-2.41063","style":"overflow:visible"})
            inkex.etree.SubElement( marker, inkex.addNS("path","svg"), 
                    {    "d":"m 6.55552,-2.41063 0,0 L 13.11104,0 c -1.0473,-1.42323 -1.04126,-3.37047 0,-4.82126",
                        "style": "fill:#000044; fill-rule:evenodd;stroke-width:0.62500000;stroke-linejoin:round;"    }
                )
                
        if group is None:
            group = self.layers[0]

        for seg in path:
            if seg[0] == "line":
                typ, start, end = seg
                inkex.etree.SubElement( group, inkex.addNS("path","svg"), {
                                                                            "d": "M %f,%f L %f,%f" % (start[0], start[1], end[0], end[1]),
                                                                            "style": line_style,
                                                                            "rightgcode": "Preview"
                                                                          })
                
            elif seg[0] == "arc":
                typ, cx, cy, r, start, end = seg
                
                start_angle = atan2(-cy, -cx) % pi2

                end_angle = atan2(end[1] - start[1] - cy, end[0] - start[0] - cx) % pi2
                
                radius = abs(r)
                
                if r < 0:
                    tmp = start_angle
                    start_angle = end_angle
                    end_angle = tmp
                
                inkex.etree.SubElement(    group, inkex.addNS('path','svg'), 
                         {
                            'style': arc_style_r if r<0 else arc_style,
                             inkex.addNS('cx','sodipodi'):        str(cx+start[0]),
                             inkex.addNS('cy','sodipodi'):        str(cy+start[1]),
                             inkex.addNS('rx','sodipodi'):        str(radius),
                             inkex.addNS('ry','sodipodi'):        str(radius),
                             inkex.addNS('start','sodipodi'):     str(start_angle),
                             inkex.addNS('end','sodipodi'):       str(end_angle),
                             inkex.addNS('open','sodipodi'):    'true',
                             inkex.addNS('type','sodipodi'):    'arc',
                            "rightgcode": "Preview",
                        })
    
    def gen_line(self, minX, minY, maxX, maxY, startX, startY, endX, endY, distance, angle, number = None):
        angle = copysign(abs(angle)%180.0, angle)
        
        full_width = abs(maxX - minX)
        full_height = abs(maxY - minY)
        
        if minX > maxX:
            tmp = minX
            minX = maxX
            maxX = tmp

        if minY > maxY:
            tmp = minY
            minY = maxY
            maxY = tmp
        
        if startX > endX:
            tmp = startX
            startX = endX
            endX = tmp

        if startY > endY:
            tmp = startY
            startY = endY
            endY = tmp

        width = endX - startX
        height = endY - startY

        sin_term = -sin(angle*pi/180)
        cos_term = cos(angle*pi/180)

        x = 0
        y = 0

        lines = []

        if abs(angle) > 45 and abs(angle) <135:
#            startX -= x_len
#            endX += x_len

            x_dist = distance/sin_term
            x_len = full_height*cos_term/sin_term

            startX = min(startX+x_len, startX-x_len)
            endX = max(endX+x_len, endX-x_len)

            if number is None:
                start_number = int(ceil((startX - minX)/x_dist))
                end_number = int(ceil((endX - minX)/x_dist))
                
                if start_number > end_number:
                    tmp = end_number
                    end_number = start_number
                    start_number = tmp
            
#                sys.stderr.write("start=%d, end=%d, startY=%f, endY=%f\n" % (start_number, end_number, startY, endY))
            
                for number in xrange(start_number, end_number+1):
                    x = minX + x_dist*number
                    y = minY - 1
            
                    x2 = x + x_len
                    y2 = maxY + 1
                    
                    lines.append({"a":[x, y], "b":[x2, y2], "kx": x2-x, "ky": y2-y, "number": number})
           
            else:
                x = minX + x_dist*number
                y = minY - 1
            
                x2 = x + x_len
                y2 = maxY + 1
                
#                sys.stderr.write("number=%d, minY=%f, maxY=%f, startY=%f, endY=%f\n" % (number, minY, maxY, startY, endY))
                
                return {"a":[x, y], "b":[x2, y2], "kx": x2-x, "ky": y2-y, "number": number}
        elif abs(angle) <= 45 or abs(angle) >=135:
#            startY -= y_len
#            endY += y_len
    
            y_dist = distance/cos_term
            y_len = -full_width*sin_term/cos_term
            
            startY = min(startY+y_len, startY-y_len)
            endY = max(endY+y_len, endY-y_len)
            
            if number is None:
                start_number = int(ceil((startY - minY)/y_dist))
                end_number = int(ceil((endY - minY)/y_dist))
#                sys.stderr.write("%f, %f\n" % (startY, endY))
#                sys.stderr.write("start=%d, end=%d, y_dist=%f\n" % (start_number, end_number, y_dist))
 
                if start_number > end_number:
                    tmp = end_number
                    end_number = start_number
                    start_number = tmp
            
                for number in xrange(start_number, end_number+1):
                    x = maxX + 1
                    y = minY + y_dist*number
            
                    x2 = minX - 1
                    y2 = y + y_len
                    
                    lines.append({"a":[x, y], "b":[x2, y2], "kx": x2-x, "ky": y2-y, "number": number})
            else:
                x = maxX + 1
                y = minY + y_dist*number
            
                x2 = minX - 1
                y2 = y + y_len
            
                return {"a":[x, y], "b":[x2, y2], "kx": x2-x, "ky": y2-y, "number": number}
        
        return lines
        
    def add_element(self, group, d):
       
            stroke_width = min(0.3, max(0.1, self.options.distance/2000))
            ele = inkex.etree.Element('{http://www.w3.org/2000/svg}path')
            group.append(ele)
            ele.set('d', d)
            ele.set("style", "fill:none;stroke:#ff0000;stroke-opacity:1;stroke-width:%fpx;" % stroke_width) #(min(0.3, self.options.distance/2000),))
        
    def effect(self):
        stroke_width = min(0.3, max(0.1, self.options.distance/2000))
        self.get_info()
		
	if self.options.add_numeric_suffix_to_filename :
            dir_list = os.listdir(self.options.directory)
            if "." in self.options.file : 
                r = re.match(r"^(.*)(\..*)$",self.options.file)
                ext = r.group(2)
                name = r.group(1)
            else:     
                ext = ""
                name = self.options.file
            max_n = 0
            for s in dir_list :
                r = re.match(r"^%s_0*(\d+)%s$"%(re.escape(name),re.escape(ext) ), s)
                if r :
                    max_n = max(max_n,int(r.group(1)))
            filename = name + "_" + ( "0"*(4-len(str(max_n+1))) + str(max_n+1) ) + ext
            self.options.file = filename
        
        if self.selected_paths == {} :
            paths=self.paths
        else :
            paths = self.selected_paths
            
        
            group = inkex.etree.SubElement( self.document.getroot(), inkex.addNS('g','svg'), {"GcodeTest": "Preview group"} )
            
        
            if self.options.paths:
                self.contours_generator = GcodeGenerator(self.options.laser_on, self.options.laser_off, self.options.laser_power, self.options.laser_on_pause, self.options.laser_off_pause, self.options.feed_speed, self.options.job_preambule, self.options.job_postambule, self.options.path_preambule, self.options.path_postambule, self.options.delete_comments)
                self.contours_generator.start(self.document.getroot().get('height'))
            
            if self.options.filling:
                self.fills_generator = GcodeGenerator(self.options.laser_on, self.options.laser_off, self.options.laser_power, self.options.laser_on_pause, self.options.laser_off_pause, self.options.feed_speed, self.options.job_preambule, self.options.job_postambule, self.options.path_preambule, self.options.path_postambule, self.options.delete_comments)
                self.fills_generator.start(self.document.getroot().get('height'))
        
        objects = []
        ids = []
        minX = float("inf")
        maxX = -float("inf")
        minY = float("inf")
        maxY = -float("inf")
        
        for layer in self.layers :
            if layer in paths :
                for path in paths[layer] :
                    d = path.get("d")
                    
                    csp = cubicsuperpath.parsePath(d)
                    csp = self.apply_transforms(path, csp)
                    
                    segments = []
                    
                    for seg in csp:
                        approx = []
                        beziers = []
                        prev_point = []
                        needs_correction = False
                        correction_point = []
                        for i in xrange(len(seg)-1):
                            point = [seg[i][1], seg[i][2], seg[i+1][0], seg[i+1][1]]
                            if len(prev_point) > 0:
                                if point[0] == prev_point[0] and point[1] == prev_point[1] and point[2] == prev_point[2] and point[3] == prev_point[3]:
                                    continue
                            
                            if needs_correction:
                                needs_correction = False
                                point[0] = correction_point[:]
                                
                            bezier = Bezier2d(point, self.options.linearity_tolerance, self.options.length_threshold/1000)
                            
                            if len(bezier.pts) == 0:
                                needs_correction = True
                                if len(beziers) == 0:
                                    correction_point = point[0][:]
                                else:
                                    correction_point = [(point[0][0]+point[-1][0])/2, (point[0][1]+point[-1][1])/2]
                                    pts = beziers[-1].pts
                                    pts[-1] = correction_point[:]
                                    beziers[-1] = Bezier2d(pts, self.options.linearity_tolerance, self.options.length_threshold/1000)
                            else:
                                prev_point = point[:]
                                beziers.append(bezier)
                                
                            if self.options.paths:
                                approx += bezier.approximate(self.options.approximation_tolerance/1000, self.options.linearity_tolerance)
                                
                        if self.options.paths:
                                     self.draw_path(approx, stroke_width, group)
							
	                             if self.options.paths:
                                      self.contours_generator.add_path(approx)      
                        
                            
                                
                        segments.append(beziers)
        
                    id = path.get("id")
                    bbox = self.get_bbox([id])
                    objects.append({"id": id, "segments":segments, "bbox": bbox})
                    ids.append(id)
                    minX = min(minX, bbox["left"])
                    maxX = max(maxX, bbox["right"])
                    minY = min(minY, bbox["bottom"])
                    maxY = max(maxY, bbox["top"])
                    
        if (self.options.filling) and (len(ids)>0):
#            lines = self.gen_lines(minX, maxX, maxY, minY, self.options.distance, self.options.angle)
            lines = {}
            cross_lines = {}
            
            if self.options.rule == "odd":
                shift = 0
            else:
                shift = 1
            
            for obj in objects:
                for p in obj["segments"]:
                    for seg in p:
                        for line in self.gen_line(minX, minY, maxX, maxY, seg.bbox["left"], seg.bbox["bottom"], seg.bbox["right"], seg.bbox["top"], mm_to_px(self.options.distance/1000), self.options.angle):
                            ts = seg.intersect_to_line(line["a"], line["b"])

                            number = line["number"]
                    
                            if number not in lines:
                                lines[number] = []

                            #append_different(lines[number], ts)
                            lines[number]+=ts[:]
                            
                        if self.options.fill_method == "grid":
                            for line in self.gen_line(minX, minY, maxX, maxY, seg.bbox["left"], seg.bbox["bottom"], seg.bbox["right"], seg.bbox["top"], mm_to_px(self.options.distance/1000), (self.options.angle+90)%180):
                                ts = seg.intersect_to_line(line["a"], line["b"])

                                number = line["number"]
                    
                                if number not in cross_lines:
                                    cross_lines[number] = []

                                #append_different(lines[number], ts)
                                cross_lines[number]+=ts[:]
            
            d = ""
            reverse = False
            for number in sorted(lines):
                line = self.gen_line(minX, minY, maxX, maxY, minX, minY, maxX, maxY, mm_to_px(self.options.distance/1000), self.options.angle, number)
                ts = lines[number]
#                ts = [0, 1]
#                print_("number=%d/%d, ts=%s" % (number, len(lines), str(ts)))
                if len(ts) > 1:
#                    sys.stderr.write(str(len(ts)))
                    ts = sorted(ts)
                    
                    if len(ts)%2 == 1:
                        ts = ts[1:]
                    
                    if reverse:
                        ts.reverse()
                        
                    reverse = not reverse
                    
                    cnt = 0
                    while cnt < len(ts)-1-shift:
    
                        t1 = ts[cnt+shift]
                        t2 = ts[cnt+shift+1]
                        
                        line_start = (lerp(line["a"][0], line["b"][0], t1), lerp(line["a"][1], line["b"][1], t1))
                        line_end = (lerp(line["a"][0], line["b"][0], t2), lerp(line["a"][1], line["b"][1], t2))
                        
                        d += "M %f,%f " % line_start
                        d += "%f,%f " % line_end
                        
                       
                        self.fills_generator.add_path( [["line", line_start, line_end]] )
                        
                        cnt += 2
#                else:
#                    print_("Empty.")

            
                self.add_element(group, d)
                    
            d = ""
                
    #            sys.stderr.write("%s\n" % str(self.options.fill_method))
            if self.options.fill_method == "grid":
                d = ""
                reverse = False

                for number in sorted(cross_lines):
                    line = self.gen_line(minX, minY, maxX, maxY, minX, minY, maxX, maxY, mm_to_px(self.options.distance/1000), (self.options.angle+90)%180, number)
                    ts = cross_lines[number]
#                   ts = [0, 1]
#                    print_("number=%d/%d, ts=%s" % (number, len(lines), str(ts)))
                    if len(ts) > 0:
#                        sys.stderr.write(str(len(ts)))
                        ts = sorted(ts)
                    
                        if reverse:
                            ts.reverse()
                        
                        reverse = not reverse
                    
                        cnt = 0

                        while cnt < len(ts)-1-shift:
    
                            t1 = ts[cnt+shift]
                            t2 = ts[cnt+shift+1]
                        
                            line_start = (lerp(line["a"][0], line["b"][0], t1), lerp(line["a"][1], line["b"][1], t1))
                            line_end = (lerp(line["a"][0], line["b"][0], t2), lerp(line["a"][1], line["b"][1], t2))
                            
                            d += "M %f,%f " % line_start
                            d += "%f,%f " % line_end
                            
                           
                            self.fills_generator.add_path([["line", line_start, line_end]])
                        
                            cnt += 2
#                    else:
#                        print_("Empty.")

                
                self.add_element(group, d)
                    
                d = ""

#            print_("Lines drawn.")
#        else:
#            sys.stderr.write("No IDs collected!")
        
#        if hasattr(self, "generator"):
#            f = open("/home/dingo/gcodes/5/hatching_test.gcode", "w")
#            f.write(self.generator.end())
#            f.close
        
#        print_("All done.")
    
        
            gcode = ""
            
            if (self.options.paths and self.options.filling):
                preambule = self.contours_generator.get_preambule()
                contours = self.contours_generator.get_gcode()
                fills = self.fills_generator.get_gcode()
                postambule = self.fills_generator.get_postambule()
		if self.options.delete_comments:
			gcode = "%s\n%s\n%s\n%s" % (preambule, contours, fills, postambule)
		else:
			gcode = "%s\n;(Start contours)\n%s\n;(Start fills)\n%s\n;(End fills)\n%s" % (preambule, contours, fills, postambule)
	if (self.options.filling and (not self.options.paths)):
		fills = self.fills_generator.get_gcode()
		postambule = self.fills_generator.get_postambule()
		gcode = self.fills_generator.end()
		
	if ((not self.options.filling) and self.options.paths):
		preambule = self.contours_generator.get_preambule()
		contours = self.contours_generator.get_gcode()
		gcode = self.contours_generator.end()       
	self.export_gcode(gcode)
    
e = RGcodeEffect()
e.affect()

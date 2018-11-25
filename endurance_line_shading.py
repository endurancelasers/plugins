#Version controll: last edited by 01.03.2018 8:20
#!/usr/bin/env python 
'''
Created by Danylo Horbatenko 2018, dnkxyz@gmail.com
Copyright (C) 2018 George Fomitchev, gf@endurancerobots.com

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'''

import os
import tempfile
import shutil
import subprocess
import math
import inkex
import png
import simplepath
import simplestyle

def saw(x):
  #The function returns a symmetric triangle wave with period 4 and varying between -1 and 1
  x = math.fmod(x, 4.0)
  x = math.fabs(x)
  if x > 2.0:
    y = 3 - x
  else:
    y = x - 1    
  return y
 
def square(x):
  #The function returns a square wave with period 4 and varying between -1 and 1
  x = math.fmod(x, 4.0)
  if 1.0 < x < 3.0:
    y = 1.0
  else:
    y = -1.0    
  return y

class LineShading(inkex.Effect):
    def __init__(self):
        inkex.Effect.__init__(self)
        self.OptionParser.add_option('-p', '--palette', action='store', type='string', dest='palette', help='Choose the colors...')
        self.OptionParser.add_option("--waveform",
          action="store", type="string",
          dest="waveform",
          help="Select the shape of the curve")        
        self.OptionParser.add_option("--num_lines",
          action="store", type="int",
          dest="num_lines",
          help="Number of lines")                          
        self.OptionParser.add_option("--min_period",
          action="store", type="float",
          dest="min_period",
          help="Minimum period (corresponds to black pixels)")    
        self.OptionParser.add_option("--max_period",
          action="store", type="float",
          dest="max_period",
          help="Maximum period (corresponds to white pixels)")  
        self.OptionParser.add_option("--min_amplitude",
          action="store", type="float",
          dest="min_amplitude",
          help="Minimum amplitude (corresponds to white pixels)")     
        self.OptionParser.add_option("--max_amplitude",
          action="store", type="float",
          dest="max_amplitude",
          help="Maximum amplitude (corresponds to black pixels)")   
        self.OptionParser.add_option("--gamma",
          action="store", type="float",
          dest="gamma",
          help="Maximum amplitude (corresponds to black pixels)")                                
        self.OptionParser.add_option("--line_width",
          action="store", type="float",
          dest="line_width",
          help="Line width")
        self.OptionParser.add_option("--units",
          action="store", type="string",
          dest="units",
          help="Units for line thickness")            
        self.OptionParser.add_option("--remove",
          action="store", type="inkbool",
          dest="remove",
          help="If True, source image is removed")                  
        self.OptionParser.add_option("--active-tab",
          action="store", type="string",
          dest="active-tab",
          help="The selected UI-tab when OK was pressed")

    def drawfunction(self, image_w, image_h, file):    
      reader = png.Reader(file)
      w, h, pixels, metadata = reader.read_flat()          
      matrice = [[1.0 for i in range(w)]for j in range(h)] 
      gumma = self.options.gamma
      if metadata['alpha']: 
        n = 4 
      else: 
        n = 3
      #RGB convert to grayscale 0.21R + 0.72G + 0.07B    
      for y in range(h):    
        for x in range(w):
           pixel_pos = (x + y * w)*n
           p = 1.0 - (pixels[pixel_pos]*0.21 + pixels[(pixel_pos+1)]*0.72 + pixels[(pixel_pos+2)]*0.07)/255.0
           matrice[y][x] = math.pow(p, 1.0/gumma)
      
      points = []
      step_y = image_h/h
      step_x = image_w/(w-1)
      min_amplitude = self.options.min_amplitude*step_y/2
      max_amplitude = self.options.max_amplitude*step_y/2
      min_period = self.options.min_period*step_y
      max_period = self.options.max_period*step_y
      min_frequency = 1.0/max_period   
      max_frequency = 1.0/min_period    
              
      #Sinusoidal wave (optimized)      
      if self.options.waveform == 'sin':
        for y in range(h):
          pi = math.pi
          phase = 0.0
          coord_x = 0.0
          amplitude = 0.0
          n_step = 0
          x0 = 0.0
          y0 = math.sin(phase)*(min_amplitude + (max_amplitude - min_amplitude)*matrice[y][x]) + (y+0.5)*step_y
          points.append([' M ',[x0, y0]])
          for x in range(w):
            period = min_period + (max_period - min_period)*(1-matrice[y][x])
            #period = 1.0/(min_frequency + (max_frequency - min_frequency)*(matrice[y][x]))
            d_phase = 2.0*pi/period*step_x            
            #calculate y
            if phase > 2.0*pi:
              if n_step > 0:
                x3 = coord_x
                y3 = -amplitude/n_step + (y+0.5)*step_y
                x2 = x3 - (x3-x0)*0.32
                y2 = y3
                x1 = x0 + (x3-x0)*0.34
                y1 = y0
                x0 = x3
                y0 = y3
                points.append([' C ',[x1, y1, x2, y2, x3, y3]])
                n_step = 0
                amplitude = 0
            elif phase < pi < (phase + d_phase):
              if n_step > 0:
                x3 = coord_x
                y3 = amplitude/n_step + (y+0.5)*step_y
                x2 = x3 - (x3-x0)*0.34
                y2 = y3
                x1 = x0 + (x3-x0)*0.32
                y1 = y0
                x0 = x3
                y0 = y3
                points.append([' C ',[x1, y1, x2, y2, x3, y3]])
                n_step = 0
                amplitude = 0
            phase = math.fmod(phase, 2.0*pi)            
            #calculate x
            if phase < 0.5*pi < (phase + d_phase):
              coord_x = (x - (phase - 0.5*pi)/d_phase)*step_x
            elif phase < 1.5*pi < (phase + d_phase):
              coord_x = (x - (phase - 1.5*pi)/d_phase)*step_x               
            phase += d_phase
            amplitude += (min_amplitude + (max_amplitude - min_amplitude)*matrice[y][x])
            n_step += 1         
          #add last point          
          if n_step > 0:
            phase = math.fmod(phase, 2.0*pi)
            if (0 < phase < 0.5*pi) or (pi < phase < 1.5*pi):
              x3 = (w-1)*step_x
              y3 = amplitude*math.sin(phase)/n_step + (y+0.5)*step_y
              x2 = x3
              y2 = y3
              x1 = x0 + (x3-x0)*0.33
              y1 = y0
              points.append([' C ',[x1, y1, x2, y2, x3, y3]])
            else:
              if coord_x > (w-1)*step_x:
                coord_x = (w-1)*step_x
              x3 = coord_x
              y3 = math.copysign( amplitude , math.sin(phase))/n_step + (y+0.5)*step_y
              x2 = x3 - (x3-x0)*0.32
              y2 = y3
              x1 = x0 + (x3-x0)*0.34
              y1 = y0
              points.append([' C ',[x1, y1, x2, y2, x3, y3]])
              if coord_x < (w-1)*step_x:
                x0 = x3
                y0 = y3
                x3 = (w-1)*step_x
                y3 = amplitude*math.sin(phase)/n_step + (y+0.5)*step_y
                x2 = x3
                y2 = y3
                x1 = x0 + (x3-x0)*0.33
                y1 = y0
                points.append([' C ',[x1, y1, x2, y2, x3, y3]])
            
      #Sinusoidal wave (Brute-force)      
      elif self.options.waveform == 'sin_b': 
        pi2 = math.pi*2.0
        for y in range(h):          
          phase = - pi2/4.0
          for x in range(w):
            period = min_period + (max_period - min_period)*(1-matrice[y][x])
            amplitude = min_amplitude + (max_amplitude - min_amplitude)*matrice[y][x]
            phase += pi2*step_x/period
            phase = math.fmod(phase, pi2)
            if x == 0:
              points.append([' M ',[x*step_x, amplitude*math.sin(phase) + (y+0.5)*step_y]])
            else:
              points.append([' L ',[x*step_x, amplitude*math.sin(phase) + (y+0.5)*step_y]])       
              
      #Saw wave            
      elif self.options.waveform == 'saw':
        for y in range(h):          
          phase = 0.0
          coord_x = 0.0
          amplitude = 0.0
          n_step = 0.0
          for x in range(w):
            period = min_period + (max_period - min_period)*(1-matrice[y][x])
            #period = 1.0/(min_frequency + (max_frequency - min_frequency)*(matrice[y][x]))
            d_phase = 4.0/period*step_x
            if phase > 4.0:
              coord_x = (x - (phase - 4.0)/d_phase)*step_x
            elif phase < 2.0 < (phase + d_phase):
              coord_x = (x - (phase - 2.0)/d_phase)*step_x   
            phase = math.fmod(phase, 4.0)
            if (phase < 1.0 < (phase + d_phase)) or (phase < 3.0 < (phase + d_phase)): 
              if n_step > 0:
                if coord_x == 0.0:
                  points.append([' M ',[coord_x, amplitude*square(phase - 1.0)/n_step + (y+0.5)*step_y]])
                else:
                  points.append([' L ',[coord_x, amplitude*square(phase - 1.0)/n_step + (y+0.5)*step_y]])
                n_step = 0
                amplitude = 0                          
            phase += d_phase
            n_step += 1.0
            amplitude += (min_amplitude + (max_amplitude - min_amplitude)*matrice[y][x])
          if n_step > 0:
            points.append([' L ',[(w-1)*step_x, amplitude*saw(phase - 1.0)/n_step + (y+0.5)*step_y]])

      #Square wave      
      else:
        for y in range(h):          
          phase = 0.0
          coord_x = 0.0
          amplitude = 0.0
          n_step = 0
          for x in range(w):
            period = min_period + (max_period - min_period)*(1-matrice[y][x])
            #period = 1.0/(min_frequency + (max_frequency - min_frequency)*(matrice[y][x]))
            d_phase = 4.0/period*step_x           
            if phase > 4.0:
              coord_x = (x - (phase - 4.0)/d_phase)*step_x
            elif phase < 2.0 < (phase + d_phase):
              coord_x = (x - (phase - 2.0)/d_phase)*step_x                
            phase = math.fmod(phase, 4.0)  
            if phase < 1.0 < (phase + d_phase): 
              if n_step > 0:
                if coord_x == 0.0:
                  points.append([' M ',[coord_x, amplitude/n_step + (y+0.5)*step_y]])
                else:
                  points.append([' L ',[coord_x, -amplitude/n_step + (y+0.5)*step_y]])
                  points.append([' L ',[coord_x, amplitude/n_step + (y+0.5)*step_y]])
                n_step = 0
                amplitude = 0     
            elif phase < 3.0 < (phase + d_phase): 
              if n_step > 0:
                if coord_x == 0.0:
                  points.append([' M ',[coord_x, -amplitude/n_step + (y+0.5)*step_y]])
                else:
                  points.append([' L ',[coord_x, amplitude/n_step + (y+0.5)*step_y]])
                  points.append([' L ',[coord_x, -amplitude/n_step + (y+0.5)*step_y]])
                n_step = 0
                amplitude = 0                              
            phase += d_phase
            n_step += 1
            amplitude += (min_amplitude + (max_amplitude - min_amplitude)*matrice[y][x])
          if n_step > 0:
            if 3.0 > phase > 1.0:
              points.append([' L ',[(w-1)*step_x, amplitude/n_step + (y+0.5)*step_y]])
            else:
              points.append([' L ',[(w-1)*step_x, -amplitude/n_step + (y+0.5)*step_y]])
      return points
  
    def draw_path(self, node, file):     
      newpath = inkex.etree.Element(inkex.addNS('path','svg'))
      line_width = self.options.line_width
      units = self.options.units
      s = {'stroke': '#000000', 'fill': 'none', 'stroke-linejoin': 'round', 'stroke-linecap': 'round', 'stroke-width': str(self.unittouu(str(line_width) + units))}
      newpath.set('style', simplestyle.formatStyle(s))      
      x = node.get('x')
      y = node.get('y')
      t = 'translate('+ x +','+ y +')'
      newpath.set('transform', t)
      image_w = float(node.get('width'))
      image_h = float(node.get('height'))          
      newpath.set('d', simplepath.formatPath(self.drawfunction(image_w, image_h, file)))
      newpath.set('title', 'Line_Shading')
      node.getparent().append(newpath)
      newpath.set('x', x)
       
    def export_png(self, node, file):
      image_w = float(node.get('width'))
      image_h = float(node.get('height'))
      min_period = self.options.min_period 
      max_period = self.options.min_period        
      poinnt_per_min_period = 8.0           
      current_file = self.args[-1]
      h_png = str(self.options.num_lines)
      if min_period < max_period:
        w_png = str(round(poinnt_per_min_period*image_w*float(h_png)/min_period/image_h))
      else:
        w_png = str(round(poinnt_per_min_period*image_w*float(h_png)/max_period/image_h))
      id = node.get('id') 
      cmd = ["inkscape", current_file, "--export-png", file, "-w", w_png, "-h", h_png, "--export-background", "rgb(255, 255, 255)", "--export-background-opacity", "255", "--export-id", id]
      proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      return_code = proc.wait()
      #f = proc.stdout
      #err = proc.stderr

    def effect(self):
      image_selected_flag = False
      for id, node in self.selected.iteritems():
        if node.tag == inkex.addNS('image','svg'):
          image_selected_flag = True
          tmp_dir = tempfile.mkdtemp()
          png_temp_file = os.path.join(tmp_dir, "LineShading.png")
          self.export_png(node, png_temp_file)                    
          self.draw_path(node, png_temp_file)
          shutil.rmtree(tmp_dir)
          if self.options.remove:
            node.getparent().remove(node)
          return
      if not image_selected_flag:
        inkex.errormsg(_("Please select an image"))

# Create effect instance and apply it.
effect = LineShading()
effect.affect()

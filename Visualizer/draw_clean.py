#  drawclean.py uses VPython to draw neutrons in a simulation geometry.
#  Geometry is defined by A.T.Holley.
#  ARR 1/11/12, 2/06/12, 2012-04-09

#from __future__ import division
#from visual import *
import os
#from ucn_import_lib import *
#from ucn_render_lib import *
from Tkinter import Tk, Label, Button

class MyFirstGUI:
    def __init__(self, master):
        self.master = master
        master.title("UCN Trajectory Visualizer")
        
        self.label = Label(master, text="Please select appropriate input and output files")
        self.label.pack()
        
        self.greet_button = Button(master, text="Greet", command=self.greet)
        self.greet_button.pack()
        
        self.close_button = Button(master, text="Close", command=master.quit)
        self.close_button.pack()
    
    def greet(self):
        print("Greetings!")

root = Tk()
root.geometry("500x100+500+400")
my_gui = MyFirstGUI(root)

root.lift()
root.attributes('-topmost',True)
root.after_idle(root.attributes,'-topmost',False)

root.mainloop()




#scene.width =1200
#scene.height = 800
#scene.autocenter = true


#regions = read_regions();
#connex = read_connex();
#pieces = draw_geom(regions, connex);

#pieces[-1].color=color.blue

#hitplaces = read_simfile('noscattercleantraj.txt');
#trails = draw_simfile(hitplaces,1);

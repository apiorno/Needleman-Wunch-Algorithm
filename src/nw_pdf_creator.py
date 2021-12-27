from fpdf import FPDF

class PDF(FPDF):
    
    def lines(self):
        self.set_fill_color(32.0, 47.0, 250.0) # color for outer rectangle
        self.rect(5.0, 5.0, 200.0,287.0,'DF')
        self.set_fill_color(255, 255, 255) # color for inner rectangle
        self.rect(8.0, 8.0, 194.0,282.0,'FD')
    def titles(self):
        self.set_xy(0.0,0.0)
        self.set_font('Arial', 'B', 16)
        self.set_text_color(220, 50, 50)
        self.cell(w=210.0, h=40.0, align='C', txt="Needleman-Wunch Alg GRASP results", border=0)
    def charts(self,plt,x,y):
        self.set_xy(x,y)
        self.image(plt,  link='', type='', w=700/5, h=450/5)
    def texts(self,name):
        with open(name,'rb') as xy:
            txt=xy.read().decode('latin-1')
        self.set_xy(10.0,80.0)    
        self.set_text_color(76.0, 32.0, 250.0)
        self.set_font('Arial', '', 12)
        self.multi_cell(0,10,txt)
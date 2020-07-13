from matplotlib import pyplot as plt
import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import io
import base64
from matplotlib.figure import Figure
from IPython.display import HTML

class figure_gif(object):

    def __init__(self, *args,**kwargs):
        self.frames=[]
        self.fig=self.init_figure(*args,**kwargs)
        return 

    def init_figure(self,*args, **kwargs):
        return plt.Figure(*args, **kwargs)
    
    def add_subplot(self,*args, **kwargs):
        return self.fig.add_subplot(*args, **kwargs)
    
    def add_frame(self):
        from io import BytesIO
        canvas=FigureCanvas(self.fig)
        png_output = BytesIO()
        canvas.print_png(png_output)
        self.frames.append(png_output)

    def save_as_html(self,filename,title=""):
        html=self.html(caption=title,fps=10)
        html_template="""<html><!doctype html><head>
             <script src="https://cdnjs.cloudflare.com/ajax/libs/jquery/1.12.0/jquery.min.js"></script>
            <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/4.7.0/css/font-awesome.min.css">

            </head>
            <body>\n{content}\n</body>
            </html>"""
        with open(filename,'w') as f:
            f.write(html_template.format(content=html))
        return

    def show(self,caption='',fps=10):
        html=self.html(caption=caption,fps=fps)
        display(HTML(html))
        plt.close(self.fig)

    def html(self,caption='',fps=10):
        from base64 import b64encode
        image_list=self.frames
        html_template="""
                <div class='animation' data-rate='0' data-fps={fps}  data-frame='0' data-images='[ {img_list} ]'>
                <figure><figcaption>{caption}</figcaption><img id='frame'/>
                    <div id="buttonbar" style='padding-left: 4em;'>
                    <button title='Reset' id="restart" ><i class="fa fa-undo"></i></button>
                    <button title='Step-backwards' id="rew" ><i class="fa fa-step-backward" ></i></button>
                    <button title='Pause' id="pause" ><i class="fa fa-pause" ></i></button>
                    <button title='Play' id="play" ><i class="fa fa-play" ></i></button>
                    <button title='Step-forward' id="fastFwd"><i class="fa fa-step-forward" ></i></button>
                    &nbsp; &nbsp; <label for="replay">frames/sec: </label>
                    <input style="width: 100px; display: inline" title="Adjust replay rate" type="range" min="0.1" step="0.1" max="30" value="{fps}" class="slider" id="myRange" oninput="v.value=value">
                    <input id=v value="{fps}" style="width: 3em; display: inline">
                    </div></figure><br/>
                </div>
                """
        javascript = """ <script>   $('div.animation').each(function(){
       	var rate=1;
        var $anim=$(this);
        var $fps=$anim.find('input.slider');
        var myInterval;
        var test=(function($anim){

        	var myfunc=function(loop=true){
                    var n=$anim.data('frame');
                    var images=$anim.data('images');
                    var n_img=images.length;
                    var image_src=images[n];                    
                    $anim.find('#frame').attr('src', image_src);
                    var increment=$anim.data('rate');
                    if (increment==0) return;
                    if (increment==-1 && n==0) return;
                    increment=increment%n_img
                    n=(n+increment+n_img)%n_img;
                    $anim.data('frame',n);
                    if (n==0){
                        $anim.data('rate',0);
                        return;
                    }
                    
                    var fps=$fps.val();
                    var ms=1000/fps;    
                    if (loop) { setTimeout(myfunc,ms); }                
                }
            return myfunc
        });
        var animate=test($anim)
        var draw=function(rate){ $anim.data('rate',rate); animate(false);  };
        //var timeout=function(){setTimeout(function(){play(0); },20000)};
        var play=function(rate){ $anim.data('rate',rate); animate(true);  };
        $anim.find('#rew').on('click',function(){draw(-1);});
        $anim.find('#fastFwd').on('click',function(){draw(1);});
        $anim.find('#play').on('click',function(){play(1);});
       	$anim.find('#pause').on('click',function(){draw(0);});
        $anim.find('#restart').on('click',function(){$anim.data('frame',0);draw(0);});
        draw(0);
    });</script> 
        """
        tag='"data:image/png;base64,{img_data}"\n'
        js_image_list=""
        sep=''
        for img in image_list:
            img_data = b64encode(img.getvalue()).decode('ascii')
            js_image_list+=sep+tag.format(img_data=img_data)
            sep=','
        html=html_template.format(img_list=js_image_list,caption=caption,fps=fps)+javascript
        return html
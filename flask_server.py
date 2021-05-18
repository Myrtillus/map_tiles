import os.path
from flask import Flask, send_file   # pip install flask
import pdb

# Oletuksena se, että .../tiles hakemiston alla zoomeja vastaavat hakemistot, sen jälkeen x-suunta, jonka alla y-suunnan png fileet.
# palauttaa tyhjän läpinäkyvän png tiedoston clientille, jos jotain tiettyä tiiltä ei ole laskettu

app = Flask(__name__, static_url_path='/home/ubuntu/map_tiles/tiles')

@app.route('/tiles/<zoom>/<x>/<y>', methods=['GET', 'POST'])
def tiles(zoom, y, x):
    default = 'tiles/tyhja.png' # this is a blank tile, change to whatever you want
    filename = '/home/ubuntu/map_tiles/tiles/%s/%s/%s' % (zoom, x, y)
 
    print(filename)
    if os.path.isfile(filename):
        return send_file(filename)
    else:
        return send_file(default)

@app.route('/', methods=['GET', 'POST'])
def index():
    return app.send_static_file('index.html')
    
if __name__ == '__main__':
    app.run(debug=False, host='localhost', port=9000)

from flask import Flask, request, render_template
from pep import PeptideOperations
import pep
from Bio import SeqIO
import io, os
import re

user_count = pep.visitorCountInit()

if not os.path.exists("mysite/static/serve/images"):
    os.makedirs("mysite/static/serve/images")
if not os.path.exists("mysite/static/serve/files"):
    os.makedirs("mysite/static/serve/files")


app = Flask(__name__, template_folder='./template', static_folder='./static')
app.config["DEBUG"] = False


@app.route("/")
def home():
    global user_count
    newCount=pep.updateVisitorCount(user_count)
    if newCount != None:
        user_count=newCount
    return render_template('home.html', userCount=user_count)


@app.route("/peptool")
def peptool():
    return render_template('pepTool.html')


@app.route("/blast_tool")
def blasttool():
    return render_template('blastTool.html')


@app.route("/calculate_result", methods=['POST'])
def calculate():
    pepOps = PeptideOperations()
    seq = None
    check_list= []
    result=None

    seq = request.form["seq"]
    check_list=request.form.getlist("cal")

    result = pepOps.check(seq, check_list)
    return render_template('pepResult.html', result=result )


@app.route("/calculate_blast_result", methods=['POST'])
def calculateBlast():
    result=None
    fileContent=None
    # Read the uploaded file
    fasta_file = request.files.get('fasta_file')
    if fasta_file is not None and fasta_file.filename != '':
        contents = fasta_file.read()
        # Convert the contents to a string
        fileContent = SeqIO.read(io.StringIO(contents.decode('utf-8')), 'fasta')
        result = pep.blast(fileContent)
    return render_template('blastResult.html', result=result)


@app.route("/documentation")
def documentation():
    return render_template('help.html')


@app.route("/bulk_sequence_list", methods=['GET', 'POST'])
def bulk_list():
    sequences = []

    if request.method == 'POST':
        seq_file = request.files.get('seq_file')

        # Regex to keep only valid amino acids
        valid_pattern = re.compile('[^ACDEFGHIKLMNPQRSTVWY]')

        if seq_file and seq_file.filename.endswith('.txt'):
            content = seq_file.read().decode('utf-8')
            # Split by new lines and clean each sequence
            raw_lines = [line.strip().upper() for line in content.splitlines() if line.strip()]

            for seq in raw_lines:
                cleaned_seq = valid_pattern.sub('', seq)
                if cleaned_seq:
                    sequences.append(cleaned_seq)

    return render_template('bulkList.html', sequences=sequences)


if __name__ == '__main__':
    # Gunicorn configuration
    gunicorn_options = {
        'bind': '0.0.0.0:8000',
        'workers': 4,
        'threads': 2
    }

    # Start the Gunicorn server
    from gunicorn.app.base import BaseApplication
    class StandaloneApplication(BaseApplication):
        def __init__(self, app, options=None):
            self.options = options or {}
            self.application = app
            super(StandaloneApplication, self).__init__()

        def load_config(self):
            config = {key: value for key, value in self.options.items() if key in self.cfg.settings and value is not None}
            for key, value in config.items():
                self.cfg.set(key.lower(), value)

        def load(self):
            return self.application

    StandaloneApplication(app, gunicorn_options).run()
'''
if __name__ == '__main__':
    app.run()
'''





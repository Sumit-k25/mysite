from flask import Flask, request, render_template, jsonify, send_file
from pep import PeptideOperations
import pep

try:
    from Bio import SeqIO
except ImportError:
    SeqIO = None

import io, os, time, textwrap
import re


def _cleanup_old_images(max_age_seconds: int = 300):
    """Remove old generated plot PNGs to prevent disk buildup.

    This is used in endpoints that generate plots for temporary rendering.
    """
    image_dir = os.path.join(app.static_folder, 'serve', 'images')
    now = time.time()

    try:
        for fname in os.listdir(image_dir):
            if not fname.lower().endswith('.png'):
                continue
            path = os.path.join(image_dir, fname)
            try:
                if now - os.path.getmtime(path) > max_age_seconds:
                    os.remove(path)
            except OSError:
                pass
    except OSError:
        pass

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

    if SeqIO is None:
        # Biopython missing; render error page so user knows what to install.
        result = "Biopython is required for BLAST functionality. Install it with: pip install biopython"
        return render_template('blastResult.html', result=result)

    # Read the uploaded file
    fasta_file = request.files.get('fasta_file')
    if fasta_file is not None and fasta_file.filename != '':
        contents = fasta_file.read()
        # Convert the contents to a string
        fileContent = SeqIO.read(io.StringIO(contents.decode('utf-8')), 'fasta')
        result = pep.blast(fileContent)
    return render_template('blastResult.html', result=result)


@app.route("/calculate_result_json", methods=['POST'])
def calculate_result_json():
    """Returns the same peptide operation result as JSON.

    Used by the bulk sequence list page to show a modal and provide a download
    link without navigating away from the page.

    If an error occurs during the calculation, returns a JSON payload with an
    "error" field so the client can show a more helpful message.
    """
    try:
        pepOps = PeptideOperations()

        # Support both form-encoded and JSON payloads.
        data = request.get_json(silent=True) or request.form
        seq = data.get("seq")
        check_list = data.get("cal") or []

        # Ensure check_list is always a list when provided as a single string.
        if isinstance(check_list, str):
            check_list = [check_list]

        # Clean up old temporary images before generating new ones.
        _cleanup_old_images()

        result = pepOps.check(seq, check_list)

        # Clean up any plot images created during calculation to avoid disk bloat.
        image_dir = os.path.join(app.static_folder, 'serve', 'images')
        for entry in result:
            if isinstance(entry, dict):
                for v in entry.values():
                    if isinstance(v, str) and v.lower().endswith('.png'):
                        try:
                            os.remove(os.path.join(image_dir, v))
                        except OSError:
                            pass

        return jsonify(result)
    except Exception as e:
        # Returning 500 with JSON error helps frontend surface the underlying issue.
        return jsonify({"error": str(e)}), 500


def _format_result_text(result: list) -> str:
    if not isinstance(result, list):
        return jsonify(result)

    lines = []
    for entry in result:
        if isinstance(entry, dict):
            for k, v in entry.items():
                lines.append(f"{k}:")
                if isinstance(v, str):
                    lines.append(f"  {v}")
                elif isinstance(v, dict):
                    for subk, subv in v.items():
                        lines.append(f"  - {subk}: {subv}")
                elif isinstance(v, list):
                    for item in v:
                        lines.append(f"  - {item}")
                else:
                    lines.append(f"  {v}")
            lines.append("")
        else:
            lines.append(str(entry))

    # Wrap long lines for better PDF layout
    wrapped = []
    for line in lines:
        wrapped.extend(textwrap.wrap(line, width=90) or [""])
    return "\n".join(wrapped)


@app.route("/generate_result_pdf", methods=['POST'])
def generate_result_pdf():
    """Generate a PDF for the same result that /calculate_result_json returns."""
    pepOps = PeptideOperations()
    seq = request.form.get("seq") or (request.json and request.json.get("seq"))
    check_list = request.form.getlist("cal") if request.form else (request.json and request.json.get("cal", []))

    # Clean up any old temporary images before starting.
    _cleanup_old_images()

    result = pepOps.check(seq, check_list)
    text = _format_result_text(result)

    # Determine any images produced by the calculation so they can be cleaned up.
    image_dir = os.path.join(app.static_folder, 'serve', 'images')
    image_file_names = []
    for entry in result:
        if isinstance(entry, dict):
            for v in entry.values():
                if isinstance(v, str) and v.lower().endswith('.png'):
                    image_file_names.append(v)

    used_image_paths = [os.path.join(image_dir, name) for name in image_file_names if os.path.exists(os.path.join(image_dir, name))]

    # Create a PDF in memory using reportlab for better formatting.
    buffer = io.BytesIO()

    try:
        from reportlab.lib import colors
        from reportlab.lib.pagesizes import A4
        from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
        from reportlab.platypus import Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle

        def _format_value(val):
            if isinstance(val, dict):
                return '<br/>'.join(f'• {k}: {v}' for k, v in val.items())
            if isinstance(val, list):
                return '<br/>'.join(f'• {v}' for v in val)
            return str(val)

        # Build table rows from the result list.
        table_data = [[Paragraph('<b>Feature</b>', getSampleStyleSheet()['Heading6']), Paragraph('<b>Value</b>', getSampleStyleSheet()['Heading6'])]]
        for entry in result:
            if isinstance(entry, dict):
                for k, v in entry.items():
                    table_data.append([Paragraph(str(k), getSampleStyleSheet()['BodyText']), Paragraph(_format_value(v), getSampleStyleSheet()['BodyText'])])
            else:
                table_data.append([Paragraph(str(entry), getSampleStyleSheet()['BodyText']), ''])

        doc = SimpleDocTemplate(buffer, pagesize=A4, rightMargin=30, leftMargin=30, topMargin=30, bottomMargin=30)
        story = []

        title_style = ParagraphStyle('title', parent=getSampleStyleSheet()['Heading1'], alignment=1, spaceAfter=6)
        meta_style = ParagraphStyle('meta', parent=getSampleStyleSheet()['BodyText'], alignment=1, spaceAfter=12, fontSize=9)

        story.append(Paragraph('PepAnalyzer Results', title_style))
        story.append(Paragraph(f'Sequence: {seq}', meta_style))
        story.append(Paragraph(f'Generated: {time.strftime("%Y-%m-%d %H:%M:%S")}', meta_style))
        story.append(Spacer(1, 12))

        tbl = Table(table_data, colWidths=[150, 350])
        tbl.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#e8e8e8')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.black),
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('VALIGN', (0, 0), (-1, -1), 'TOP'),
            ('FONTNAME', (0, 0), (-1, -1), 'Helvetica'),
            ('FONTSIZE', (0, 0), (-1, -1), 9),
            ('BOTTOMPADDING', (0, 0), (-1, 0), 8),
            ('TOPPADDING', (0, 0), (-1, 0), 8),
            ('BACKGROUND', (0, 1), (-1, -1), colors.whitesmoke),
            ('LEFTPADDING', (0, 0), (-1, -1), 6),
            ('RIGHTPADDING', (0, 0), (-1, -1), 6),
            ('BOX', (0, 0), (-1, -1), 0.5, colors.grey),
            ('GRID', (0, 0), (-1, -1), 0.25, colors.grey),
        ]))

        story.append(tbl)

        # Attach any generated images (plots/figures) present in the result.
        try:
            from reportlab.platypus import Image
            for img_path in used_image_paths:
                story.append(Spacer(1, 12))
                story.append(Paragraph(os.path.basename(img_path).replace('_', ' ').rsplit('.', 1)[0].title(), getSampleStyleSheet()['Heading6']))
                story.append(Spacer(1, 6))
                story.append(Image(img_path, width=400, height=400, kind='proportional'))
        except Exception:
            # If image embedding fails, ignore and continue with PDF generation.
            pass

        doc.build(story)

        buffer.seek(0)
        filename = f"analysis_{int(time.time())}.pdf"
        return send_file(buffer, mimetype='application/pdf', as_attachment=True, download_name=filename)

    except ImportError:
        # If reportlab isn't installed, fallback to plain-text PDF in matplotlib.
        try:
            from matplotlib.backends.backend_pdf import PdfPages
            import matplotlib.pyplot as plt

            with PdfPages(buffer) as pdf:
                fig = plt.figure(figsize=(8.27, 11.69))
                fig.text(0.01, 0.99, text, va='top', family='monospace', fontsize=10)
                pdf.savefig(fig)
                plt.close(fig)

            buffer.seek(0)
            filename = f"analysis_{int(time.time())}.pdf"
            return send_file(buffer, mimetype='application/pdf', as_attachment=True, download_name=filename)
        except Exception:
            buffer = io.BytesIO(text.encode('utf-8'))
            buffer.seek(0)
            return send_file(buffer, mimetype='text/plain', as_attachment=True, download_name=f"analysis_{int(time.time())}.txt")
    except Exception:
        buffer = io.BytesIO(text.encode('utf-8'))
        buffer.seek(0)
        return send_file(buffer, mimetype='text/plain', as_attachment=True, download_name=f"analysis_{int(time.time())}.txt")
    finally:
        for p in used_image_paths:
            try:
                os.remove(p)
            except OSError:
                pass


@app.route("/documentation")
def documentation():
    return render_template('help.html')


@app.route("/bulk_sequence_list", methods=['GET', 'POST'])
def bulk_list():
    sequences = []
    selected_features = None

    if request.method == 'POST':
        seq_file = request.files.get('seq_file')
        selected_features = request.form.getlist('cal')

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

    return render_template('bulkList.html', sequences=sequences, selected_features=selected_features)


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





from flask import Flask, render_template, request
from Executors.executor import executor

app = Flask(__name__)


@app.route("/")
def home():
    return render_template("home.html")


@app.route("/about")
def about():
    return render_template("about.html")


@app.route("/leo")
def salvador():
    return "Hello, Leo. Who's a good boy?"


@app.route('/c_elegans_orthologs')
def get_c_elegans_orthologs_input():
    return render_template('query_form.html')


@app.route('/c_elegans_orthologs', methods=['POST'])
def return_c_elegans_orthologs():
    text = request.form['text']
    text = text.replace(" ", "")
    human_genes = text.split(",")
    print("human genes:", human_genes)
    genes_in_names = request.form.get('type_select')
    print("genes in names:", genes_in_names)
    try:
        results, error = executor.find_me_orthologs(human_genes=human_genes, genes_in_names=genes_in_names)
        if not results:
            query = ", ".join(human_genes)
            return render_template('failure_response.html', query=query, error=error)
    except:
        query = ", ".join(human_genes)
        error = "Something went wrong, please check your service log"
        return render_template('failure_response.html', query=query, error=error)
    return render_template('orthologs_table_response.html', results=results)


@app.route('/human_orthologs')
def get_human_orthologs_input():
    return render_template('query_form.html')


@app.route('/human_orthologs', methods=['POST'])
def return_human_orthologs():
    text = request.form['text']
    text = text.replace(" ", "")
    c_elegans_genes = text.split(",")
    print("C.elegans genes:", c_elegans_genes)
    genes_in_names = request.form.get('type_select')
    print("genes in names:", genes_in_names)
    try:
        results, error = executor.find_me_orthologs_for_worm(c_elegans_genes, genes_in_names=genes_in_names)
        if not results:  # no orthologs
            query = ", ".join(c_elegans_genes)
            return render_template('failure_response.html', query=query, error=error)
    except:
        query = ", ".join(c_elegans_genes)
        error = "Something went wrong, please check your service log"
        return render_template('failure_response.html', query=query, error=error)
    return render_template('orthologs_table_response.html', results=results)


@app.route('/variants')
def get_variants_input():
    return render_template('query_form.html')


@app.route('/variants', methods=['POST'])
def return_variants_data():
    text = request.form['text']
    text = text.replace(" ", "")
    print(text)
    genes_and_variants = executor.parse_input(text)
    print(genes_and_variants)
    try:
        results, error = executor().get_variants_data_for_server(genes_and_variants)
        if not results:
            query = executor.dictionary_output_parser(genes_and_variants)
            return render_template('failure_response.html', query=query, error=error)
    except:
        query = executor.dictionary_output_parser(genes_and_variants)
        return render_template('failure_response.html', query=query, error="Unknown")
    return render_template('variants_table_response.html', results=results)


@app.route("/post_field", methods=["GET", "POST"])
def need_input():
    for key, value in request.form.items():
        print("key: {0}, value: {1}".format(key, value))

if __name__ == "__main__":
    app.run('0.0.0.0', debug=True)

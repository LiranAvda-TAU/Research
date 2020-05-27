from flask import Flask, render_template, request
from Code.CRISPR.Enum.AminoAcid import AminoAcid
from Executors.executor import executor
from Code.CRISPR.CrisprPlanner import CrisprPlanner

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
    sources_bar = request.form.get('sources_bar')
    print("genes in names:", genes_in_names)
    try:
        results, error = executor.find_me_orthologs_for_human(human_genes=human_genes,
                                                              genes_in_names=genes_in_names,
                                                              sources_bar=sources_bar)
    except Exception as e:
        query = ", ".join(human_genes)
        error = "Something went wrong: " + str(e)
        return render_template('failure_response.html', query=query, error=error)
    true_results, false_results = results
    return render_template('orthologs_table_response.html', true_results=true_results, false_results=false_results)


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
    sources_bar = request.form.get('sources_bar')
    print("genes in names:", genes_in_names)
    try:
        results, error = executor.find_me_orthologs_for_worm(worm_genes=c_elegans_genes,
                                                             genes_in_names=genes_in_names,
                                                             sources_bar=sources_bar)
    except Exception as e:
        query = ", ".join(c_elegans_genes)
        error = "Something went wrong: " + str(e)
        return render_template('failure_response.html', query=query, error=error)
    true_results, false_result = results
    return render_template('orthologs_table_response.html', true_results=true_results, false_results=false_result)


@app.route('/variants')
def get_variants_input():
    return render_template('query_form.html')


@app.route('/variants', methods=['POST'])
def return_variants_data():
    text = request.form['text']
    text = text.replace(" ", "")
    print(text)
    genes_and_variants = executor.parse_input(text)
    sources_bar = request.form.get('sources_bar')
    print(genes_and_variants)
    try:
        true_results, false_results = executor().get_variants_data_for_server(genes_and_variants, sources_bar)
    except Exception as e:
        query = executor.dictionary_output_parser(genes_and_variants)
        return render_template('failure_response.html', query=query, error=e)
    return render_template('variants_table_response.html', true_results=true_results, false_results=false_results)


@app.route("/crispr")
def crispr_planner():
    return render_template("crispr_form.html")


@app.route('/crispr', methods=['POST'])
def return_crispr_plan():
    worm_gene_name = request.form['name']
    nt_seq = request.form['seq']
    site = int(request.form['site'])
    from_aa = AminoAcid[request.form.get('from_aa')]
    to_aa = AminoAcid[request.form.get('to_aa')]
    print("CRISPR Request:", worm_gene_name, site, nt_seq, from_aa, to_aa)
    try:
        result, error = CrisprPlanner(gene_name=worm_gene_name,
                                      aa_mutation_site=site,
                                      sense_strand=nt_seq).plan_my_crispr(from_aa=from_aa, to_aa=to_aa)
        if not result:
            return render_template('failure_response.html', query=worm_gene_name, error=error)
    except Exception as e:
        error = "Something went wrong: " + str(e)
        return render_template('failure_response.html', query=worm_gene_name, error=error)
    return render_template('crispr_response.html', result=result)


@app.route("/faq")
def faq():
    return render_template("faq.html")


@app.route("/references")
def references():
    return render_template("references.html")


if __name__ == "__main__":
    app.run(host="127.0.0.1", debug=True, port=5005)

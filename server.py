#!/usr/bin/python
from http.server import BaseHTTPRequestHandler, HTTPServer
from Executors.executor import executor

PORT_NUMBER = 8080


# This class will handles any incoming request from the browser
class HomologyHandler(BaseHTTPRequestHandler):

    def parse_input(self):
        # path = "/Variants;TCP1:{Asn284Ser,Ala453Glu},DAP3:{Leu138Phe,Glu369Lys}"
        dic = {}
        input_values = self.path[self.path.find(";")+1:]
        while input_values:
            gene = input_values[:input_values.find(":")]
            input_values = input_values[len(gene)+1:]
            variations_str = input_values[input_values.find("[")+1:input_values.find("]")]
            variations = variations_str.split(",")
            dic[gene] = variations
            input_values = input_values[len(variations_str)+3:]
        return dic

    # Handler for the GET requests
    def do_GET(self):
        self.send_response(200)
        self.send_header('Content-type', 'text/html')
        self.end_headers()

        if self.path == "/favicon.ico":
            return
        program = self.path[1:self.path.find(";")]
        print(program)
        if program == "Orthologs":
            human_genes = self.path.split(";")[1:]
            print(str(human_genes))
            try:
                c_elegans_homologs = executor.find_me_orthologs(human_genes)
            except SystemExit:
                c_elegans_homologs = "Couldn't find orthologs for " + str(human_genes)
            self.wfile.write(str(c_elegans_homologs).encode())
        elif program == "Variants":
            genes_and_variants = self.parse_input()
            try:
                results = executor.get_variants_data_for_server(genes_and_variants)
            except SystemExit:
                results = "Couldn't find data for " + str(genes_and_variants)
            self.wfile.write(results.encode())
        else:
            self.wfile.write("Program name is invalid, can only accept \"Orthologs\" or \"Variants\". Please try again".encode())
        return



try:
    # Create a web server and define the handler to manage the
    # incoming request
    server = HTTPServer(('', PORT_NUMBER), HomologyHandler)
    print('Started httpserver on port ', PORT_NUMBER)

    # Wait forever for incoming htto requests
    server.serve_forever()

except KeyboardInterrupt:
    print('^C received, shutting down the web server')
    server.socket.close()

#!/usr/bin/python
from http.server import BaseHTTPRequestHandler, HTTPServer
from Executors.executor import executor

PORT_NUMBER = 8080


# This class will handles any incoming request from the browser
class HomologyHandler(BaseHTTPRequestHandler):

    # Handler for the GET requests
    def do_GET(self):
        self.send_response(200)
        self.send_header('Content-type', 'text/html')
        self.end_headers()
        # Send the html message
        if self.path == "/favicon.ico":
            return
        human_genes = self.path[1:].split(",")
        print(str(human_genes))
        try:
            c_elegans_homologs = executor.find_me_orthologs(human_genes)
        except SystemExit:
            c_elegans_homologs = "No found orthologs for " + human_genes
        self.wfile.write(str(c_elegans_homologs).encode())
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

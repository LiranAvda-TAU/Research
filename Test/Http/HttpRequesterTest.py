from Code.Http.HttpRequester import HttpRequester


class HttpRequesterTest:

    @staticmethod
    def makeRequestTest():
        url = "http://ortholist.shaye-lab.org/"
        geneName = "usp11"
        httpRequester = HttpRequester(url, geneName)
        httpRequester.make_request()

    @staticmethod
    def test_get_transcript():
        print("###TEST: test_get_transcript###")
        list_of_gene_ids = 'WBGene00000018', 'WBGene00018050', 'WBGene00018064', 'WBGene00018072', 'WBGene00018119',\
                           'WBGene00200942', 'WBGene00200958', 'WBGene00200968', 'WBGene00201208', 'WBGene00198430',\
                           'WBGene00271642', 'WBGene00019039', 'WBGene00018933', 'WBGene00000377'
        for gene_id in list_of_gene_ids:
            print("for gene id:", gene_id)
            HttpRequester.get_transcript(gene_id)


# TESTING #

# HttpRequesterTest.makeRequestTest()
HttpRequesterTest.test_get_transcript()

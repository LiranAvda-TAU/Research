from Code.Http.HttpRequester import HttpRequester


class HttpRequesterTest:

    @staticmethod
    def makeRequestTest():
        url = "http://ortholist.shaye-lab.org/"
        geneName = "usp11"
        httpRequester = HttpRequester(url, geneName)
        httpRequester.make_request()




# TESTING #

HttpRequesterTest().makeRequestTest()

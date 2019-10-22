from django.shortcuts import render
from django.http import HttpResponse
from rest_framework import viewsets
from .serializers import GeneSerializer
from .models import Gene
# Create your views here.
# def hi(request):
#     return HttpResponse('<h1>THIS IS MY HOMEPAGE</h1>')


class GeneViewSet(viewsets.ModelViewSet):
    queryset = Gene.objects.all().order_by('name')
    serializer_class = GeneSerializer

homology_results = [
    {
        'name': 'TCP1',
        'orthologs': 'cct-1'
    },
    {
        'name': 'NEO1',
        'orthologs': 'unc-40'
    }

]
def home(request):
    results = {
        'homology_results': homology_results
    }
    return render(request, 'RESEARCH_APP/home-page.html', results)

def about(request):
    return render(request, 'RESEARCH_APP/about.html')

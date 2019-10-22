from rest_framework import serializers
from .models import Gene

class GeneSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = Gene
        fields = ('name', 'orthologs')
from django.db import models

# Create your models here.

class Gene(models.Model):
    name = models.CharField(max_length=60)
    orthologs = models.CharField(max_length=60)

    def __str__(self):
        return self.name
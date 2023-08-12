from django.db import models

# Create your models here.
class Drug(models.Model):
    name = models.CharField(max_length=255)
    chembl_id = models.CharField(max_length=20, unique=True)
    smiles = models.TextField()

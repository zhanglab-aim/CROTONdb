from django.contrib import admin
from .models import Gene, GeneInterval, CrossRef, CrossRefDB



class GeneAdmin(admin.ModelAdmin):    
    list_display = ('entrez', 'standard_name', 'ensemble_str')
    

class GeneIntervalAdmin(admin.ModelAdmin):    
    list_display = ('chromosome', 'gene', 'gene_entrez','start','end','difference')
    list_filter = ('chromosome',)

admin.site.register(Gene,GeneAdmin)
admin.site.register(GeneInterval,GeneIntervalAdmin)
admin.site.register(CrossRef)
admin.site.register(CrossRefDB)
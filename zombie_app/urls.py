from django.urls import path
from zombie_app.views import *
from . import views
from django.conf import settings\

urlpatterns = [
    path('', views.home , name='home'),
    path('download-pdf/', views.download_pdf, name='download_pdf'),
    # path('result/', views.result, name='result'),
    # path('view-pdf/', views.view_pdf, name='view_pdf'),
]
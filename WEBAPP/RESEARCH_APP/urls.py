from django.urls import path, include
from . import views
from rest_framework import routers

router = routers.DefaultRouter()
router.register(r'genes', views.GeneViewSet)

urlpatterns = [
    path('', views.home, name="home-page"),
    path('about/', views.about, name="about"),
    path('router/', include(router.urls)),
    path('api-auto/', include('rest_framework.urls', namespace='rest_framework'))
]

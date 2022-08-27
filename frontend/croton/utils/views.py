from django.views.generic.detail import BaseDetailView, SingleObjectTemplateResponseMixin
from django.views.generic.list import BaseListView, MultipleObjectTemplateResponseMixin

from utils.mixins import JSONResponseMixin

class HybridDetailView(JSONResponseMixin, SingleObjectTemplateResponseMixin, BaseDetailView):
    def render_to_response(self, context):
        if self.request.is_ajax():
            obj = context['object']
            if callable(getattr(obj, 'as_dict', None)):
                obj = obj.as_dict()
            return JSONResponseMixin.render_to_response(self, obj)
        else:
            return SingleObjectTemplateResponseMixin.render_to_response(self, context)

class HybridListView(JSONResponseMixin, MultipleObjectTemplateResponseMixin, BaseListView):
    def render_to_response(self, context):
        if self.request.is_ajax():
            obj_list = context['object_list']
            obj = [obj.as_dict() if callable(getattr(obj, 'as_dict', None)) else obj for obj in obj_list]
            return JSONResponseMixin.render_to_response(self, obj)
        else:
            return MultipleObjectTemplateResponseMixin.render_to_response(self, context)


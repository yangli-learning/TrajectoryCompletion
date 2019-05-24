#include "common.h"
#include "renderable.h"
#include "osg_utility.h"

#include "pick_handler.h"

PickHandler::PickHandler(void) {
}

PickHandler::~PickHandler(void) {
}

void PickHandler::getUsage(osg::ApplicationUsage &usage) const {
  usage.addKeyboardMouseBinding("Click+Modifier Key", "Pick Renderable Object");
  return;
}

bool PickHandler::handle(const osgGA::GUIEventAdapter& ea,osgGA::GUIActionAdapter& aa) {
  switch(ea.getEventType()) {
  case(osgGA::GUIEventAdapter::PUSH): {
    osgViewer::View* view = dynamic_cast<osgViewer::View*>(&aa);
    if (view == nullptr)
      return false;

    int mod_key_mask = ea.getModKeyMask();

	const bool ctrlKey((mod_key_mask   & osgGA::GUIEventAdapter::MODKEY_CTRL) != 0 );
    const bool shiftKey((mod_key_mask  & osgGA::GUIEventAdapter::MODKEY_SHIFT) != 0 );
    const bool altKey((mod_key_mask  & osgGA::GUIEventAdapter::MODKEY_ALT) != 0 );


    PickMode pick_mode = PickMode::CTRLSHIFT;

	if (ctrlKey && shiftKey)
		pick_mode = PickMode::CTRLSHIFT;
	else if (altKey && shiftKey)
		pick_mode = PickMode::ALTSHIFT;
    else if (ctrlKey)
      pick_mode = PickMode::CTRL;
    else if (shiftKey)
      pick_mode = PickMode::SHIFT;
    else if (altKey)
      pick_mode = PickMode::ALT;
    else
		return false;

    osgUtil::LineSegmentIntersector::Intersection intersection;
    osg::NodePath node_path;

    Renderable* renderable = OSGUtility::computeIntersection<Renderable>(view, ea, intersection, node_path);
    //Renderable* renderable = OSGUtility::computePointIntersection<Renderable>(view, ea, intersection, node_path);
    if (renderable == nullptr)
      return false;
	osg::Node *lastNode;
	if (!node_path.empty()){
		lastNode = node_path.back();
	}
	osg::Vec3d localPt = intersection.getLocalIntersectPoint();
    renderable->pickEvent(pick_mode,localPt[0],localPt[1]);
    return true;
  }
  break;
  default:
    return false;
  }

  return false;
}

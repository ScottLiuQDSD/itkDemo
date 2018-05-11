#include "SpatialObjects\include\itkSpatialObject.h"
#include "SpatialObjects\include\itkSpatialObjectTreeContainer.h"
#include "SpatialObjects\include\itkGroupSpatialObject.h"
#include "Common\include\itkLevelOrderTreeIterator.h"

void SpatialObjectHierarchy5_2()
{
	typedef itk::SpatialObject<3> SpatialObjectType;
	SpatialObjectType::Pointer object1 = SpatialObjectType::New();
	object1->GetProperty()->SetName("First Object");
	SpatialObjectType::Pointer object2 = SpatialObjectType::New();
	object2->GetProperty()->SetName("Second Object");

	object1->AddSpatialObject(object2);

	if (object2->HasParent()) {
		std::cout << "Name of parent of the object2 is :";
		std::cout << object2->GetParent()->GetProperty()->GetName() << std::endl;
	}

	SpatialObjectType::ChildrenListType *childrenList = object1->GetChildren();
	std::cout << "object1 has " << childrenList->size() << " child." << std::endl;
	SpatialObjectType::ChildrenListType::const_iterator it = childrenList->begin();
	while (childrenList->end() != it) {
		std::cout << "Name of the child of the object 1 :";
		std::cout << (*it)->GetProperty()->GetName() << std::endl;
		++it;
	}
	delete childrenList;

	object1->RemoveSpatialObject(object2);
	std::cout << "Number of children of the object1 is :";
	std::cout << object1->GetNumberOfChildren() << std::endl;
	object1->Clear();


}

void SpatialObjectTreeContainer5_3()
{
	typedef itk::GroupSpatialObject<2> NodeType;
	typedef itk::SpatialObjectTreeContainer<2> TreeType;

	NodeType::Pointer object0 = NodeType::New();
	object0->SetId(0);
	NodeType::Pointer object1 = NodeType::New();
	object1->SetId(1);
	NodeType::Pointer object2 = NodeType::New();
	object2->SetId(2);

	object0->AddSpatialObject(object1);
	object1->AddSpatialObject(object2);

	TreeType::Pointer tree = TreeType::New();
	tree->SetRoot(object0.GetPointer());
	itk::LevelOrderTreeIterator<TreeType> levelIt(tree, 10);
	levelIt.GoToBegin();
	while (!levelIt.IsAtEnd()) {
		std::cout << levelIt.Get()->GetId() << "(" << levelIt.GetLevel() << ")" << std::endl;
		++levelIt;
	}
	NodeType::Pointer object4 = NodeType::New();
	itk::PreOrderTreeIterator<TreeType> preIt(tree);
	preIt.Add(object4.GetPointer());

}

void SpatialObjectTransforms5_4()
{
	typedef itk::SpatialObject<2> SpatialObjectType;
	typedef SpatialObjectType::TransformType  TransformType;
	SpatialObjectType::Pointer object1 = SpatialObjectType::New();
	object1->GetProperty()->SetName("First Object");
	SpatialObjectType::Pointer object2 = SpatialObjectType::New();
	object2->GetProperty()->SetName("Second Object");

	object1->AddSpatialObject(object2);
	double scale[2];
	scale[0] = 2.0;
	scale[1] = 2.0;
	object2->GetIndexToObjectTransform()->SetScale(scale);

	TransformType::OffsetType Object2ToObject1Offset;
	Object2ToObject1Offset[0] = 4;
	Object2ToObject1Offset[1] = 3;
	object2->GetObjectToParentTransform()->SetOffset(Object2ToObject1Offset);
	object2->ComputeObjectToWorldTransform();

	std::cout << "object2 IndexToObject Matrix: " << std::endl;
	std::cout << object2->GetIndexToObjectTransform()->GetMatrix() << std::endl;

	std::cout << "object2 IndexToObject Offset: " << std::endl;
	std::cout << object2->GetIndexToObjectTransform()->GetOffset() << std::endl;

	std::cout << "object2 IndexToWorld Matrix: " << std::endl;
	std::cout << object2->GetIndexToWorldTransform()->GetMatrix() << std::endl;

	std::cout << "object2 IndexToWorld Offset: " << std::endl;
	std::cout << object2->GetIndexToWorldTransform()->GetOffset() << std::endl;

	TransformType::OffsetType Object1ToWorldOffset;
	Object1ToWorldOffset[0] = 3;
	Object1ToWorldOffset[1] = 3;
	object1->GetObjectToParentTransform()->SetOffset(Object1ToWorldOffset);
	object1->ComputeObjectToWorldTransform();
	std::cout << "object1 IndexToObject Matrix: " << std::endl;
	std::cout << object1->GetIndexToObjectTransform()->GetMatrix() << std::endl;

	std::cout << "object1 IndexToObject Offset: " << std::endl;
	std::cout << object1->GetIndexToObjectTransform()->GetOffset() << std::endl;

	std::cout << "object1 IndexToWorld Matrix: " << std::endl;
	std::cout << object1->GetIndexToWorldTransform()->GetMatrix() << std::endl;

	std::cout << "object1 IndexToWorld Offset: " << std::endl;
	std::cout << object1->GetIndexToWorldTransform()->GetOffset() << std::endl;

}

#include "Common\include\itkTreeContainer.h"
#include "Common\include\itkChildTreeIterator.h"
#include "Common\include\itkLeafTreeIterator.h"
#include "Common\include\itkLevelOrderTreeIterator.h"
#include "Common\include\itkInOrderTreeIterator.h"
#include "Common\include\itkPostOrderTreeIterator.h"
#include "Common\include\itkRootTreeIterator.h"
#include "Common\include\itkTreeIteratorClone.h"

LRESULT TreeContainerDemo()
{
	typedef int NodeType;
	typedef itk::TreeContainer<NodeType> TreeType;
	TreeType::Pointer tree = TreeType::New();
	tree->SetRoot(0);

	tree->Add(1, 0);
	tree->Add(2, 0);
	tree->Add(3, 0);
	tree->Add(4, 2);
	tree->Add(5, 2);
	tree->Add(6, 5);
	tree->Add(7, 1);

	itk::ChildTreeIterator<TreeType> childIt(tree);
	for (childIt.GoToBegin(); !childIt.IsAtEnd(); ++ childIt) {
		std::cout << childIt.Get() << std::endl;
	}
	childIt.GoToBegin();
	if (itk::TreeIteratorBase<TreeType>::CHILD != childIt.GetType()) {
		std::cerr << "Error: The iterator was not of type CHILD." << std::endl;
		return EXIT_FAILURE;
	}

	int oldVal = childIt.Get();
	std::cout << "The node's value is " << oldVal << std::endl;
	int newValue = 2;
	childIt.Set(newValue);
	std::cout << "Now, the node's value is " << childIt.Get() << std::endl;

	std::cout << "Is this a leaf node?" << childIt.IsLeaf() << std::endl;
	std::cout << "Is this a root node?" << childIt.IsRoot() << std::endl;
	std::cout << "Does this node have a parent?" << childIt.HasParent() << std::endl;
	std::cout << "How many children does this node have?" << childIt.CountChildren() << std::endl;
	std::cout << "Does this node have a child 1?" << childIt.HasChild(1) << std::endl;

	tree->Clear();
	itk::PreOrderTreeIterator<TreeType> it(tree);
	it.GoToBegin();
	it.Add(0);
	it.Add(1);
	it.Add(2);
	it.Add(3);
	it.GoToChild(2);
	it.Add(4);
	it.Add(5);

	itk::TreeIteratorBase<TreeType>* childItClone = childIt.Clone();
	delete childItClone;

	typedef itk::TreeIteratorBase<TreeType> IteratorType;
	typedef itk::TreeIteratorClone<IteratorType> IteratorCloneType;
	IteratorCloneType anotherChildItClone = childIt;
	for (childIt.GoToBegin(); !childIt.IsAtEnd(); ++childIt) {
		std::cout << childIt.Get();
	}
	std::cout << std::endl;

	itk::LeafTreeIterator<TreeType> leafIt(tree);
	for (leafIt.GoToBegin(); !leafIt.IsAtEnd(); ++ leafIt) {
		std::cout << leafIt.Get() << std::endl;
	}
	std::cout << std::endl;
	
	itk::LevelOrderTreeIterator<TreeType> levelIt(tree, 10, tree->GetNode(0));
	for (levelIt.GoToBegin(); !levelIt.IsAtEnd(); ++levelIt) {
		std::cout << levelIt.Get() 
			<< "(" << levelIt.GetLevel() << ") "
			<< std::endl;
	}
	std::cout << std::endl;

	itk::InOrderTreeIterator<TreeType> inOrderIt(tree);
	for (inOrderIt.GoToBegin(); !inOrderIt.IsAtEnd(); ++inOrderIt) {
		std::cout << inOrderIt.Get()
			<< std::endl;
	}
	std::cout << std::endl;
	itk::PreOrderTreeIterator<TreeType> preOrderIt(tree);
	for (preOrderIt.GoToBegin(); !preOrderIt.IsAtEnd(); ++preOrderIt) {
		std::cout << preOrderIt.Get()
			<< std::endl;
	}
	std::cout << std::endl;

	itk::PostOrderTreeIterator<TreeType> postOrderIt(tree);
	for (postOrderIt.GoToBegin(); !postOrderIt.IsAtEnd(); ++postOrderIt) {
		std::cout << postOrderIt.Get()
			<< std::endl;
	}
	std::cout << std::endl;

	itk::RootTreeIterator<TreeType> rootIt(tree, tree->GetNode(4));
	for (rootIt.GoToBegin(); !rootIt.IsAtEnd(); ++rootIt) {
		std::cout << rootIt.Get()
			<< std::endl;
	}
	std::cout << std::endl;

	return 1;
}
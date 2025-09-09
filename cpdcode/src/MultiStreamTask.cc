#include "electron_flux/MultiStreamTask.h"
#include "SniperKernel/SniperLog.h"
#include "SniperKernel/DeclareDLE.h"
#include "SniperKernel/SvcBase.h"
#include "PodioSvc/PodioInputSvc.hh"
#include <boost/python.hpp>

SNIPER_DECLARE_DLE(MultiStreamTask);

MultiStreamTask::MultiStreamTask(const std::string &name)
    : Task(name),
      m_tasks(name, "subtasks")
{
    m_tag = "MultiStreamTask";
    m_targets.push_back(&m_tasks);
}

MultiStreamTask::~MultiStreamTask()
{
}

Task *MultiStreamTask::createTask(const std::string &taskName)
{
    Task* task = new Task(taskName);
    m_tasks.append(task,true);
    task->setParent(this);
    return task;
}

Task *MultiStreamTask::addTask(Task *task)
{
    if (m_tasks.append(task, false))
    {
        task->setParent(this);
        return task;
    }
    return nullptr;
}

DLElement *MultiStreamTask::create(const std::string &type, const std::string &name)
{
    if (type == "task")
    {
        return createTask(name);
    }
    else
    {
        return Task::create(type, name);
    }
}

DLElement *MultiStreamTask::find(const std::string &name)
{
    std::string::size_type pseg = name.find(":");
    if (pseg == std::string::npos)
    {
        return Task::find(name);
    }

    Task *obj = dynamic_cast<Task *>(m_tasks.find(name.substr(0, pseg)));
    if (obj != nullptr)
    {
        return obj->find(name.substr(pseg + 1, std::string::npos));
    }

    LogWarn << "Cann't find Object " << name << std::endl;
    return nullptr;
}

void MultiStreamTask::remove(const std::string &name)
{
    std::string::size_type pseg = name.find(":");
    if (pseg == std::string::npos)
    {
        return Task::remove(name);
    }

    Task *obj = dynamic_cast<Task *>(m_tasks.find(name.substr(0, pseg)));
    if (obj != nullptr)
    {
        return obj->remove(name.substr(pseg + 1, std::string::npos));
    }

    LogWarn << "Cannot remove, no such element " << name << std::endl;
    return;
}

bool MultiStreamTask::config()
{
    bool stat = Task::config();

    for (auto obj : m_tasks.list())
    {
        Task *task = dynamic_cast<Task *>(obj);
        if (!task->Snoopy().config())
            stat = false;
    }

    if (!stat)
        m_snoopy.setErr();

    return stat;
}

bool MultiStreamTask::initialize()
{
    bool stat = true;

    if (!Task::initialize())
        stat = false;

    for (auto obj : m_tasks.list())
    {
        Task *task = dynamic_cast<Task *>(obj);
        if (!task->Snoopy().initialize())
            stat = false;
    }

    if (!stat)
        m_snoopy.setErr();

    return stat;
}

bool MultiStreamTask::finalize()
{
    bool stat = true;

    auto &tasks = m_tasks.list();
    for (auto it = tasks.rbegin(); it != tasks.rend(); ++it)
    {
        Task *task = dynamic_cast<Task *>(*it);
        if (!task->Snoopy().finalize())
            stat = false;
    }

    if (!Task::finalize())
        stat = false;

    if (!stat)
        m_snoopy.setErr();

    return stat;
}

bool MultiStreamTask::addInputFile(const std::string& fileName, const std::string& taskName)
{
    // Create a sub-task that holds this input stream

    if (fileName.empty()) {
        LogFatal << "Empty file name is not allowed. Failed to create input stream" << std::endl;
        return false;
    }

    // Create a sub-task for each input stream.
    // If task name is not provided, use fileName as the task name
    std::string streamName(fileName);
    if (!taskName.empty())  streamName = taskName;
    Task* task = this->createTask(streamName);
    task->setLogLevel(m_logLevel);
    m_streams.push_back(streamName);

    // Create PodioDataSvc and PodioInputSvc for each sub-task
    auto pds = task->createSvc("PodioDataSvc");
    auto pis = dynamic_cast<PodioInputSvc*>(task->createSvc("PodioInputSvc/InputSvc"));
    if (nullptr == pds or nullptr == pis) {
        LogFatal << "Failed to configure DataManagement services for sub-task" << std::endl;
        return false;
    }
    pis->setInputFile(fileName);
    return true;
}

bool MultiStreamTask::execute()
{
    bool stat = Task::execute();
    if (!stat) return stat;
    // After execution of TopTask
    int ret = 0;
    for (auto task: m_streams) {
        int this_ret = 0;
        try {
            this_ret = Incident::fire(*this, task); 
        }
        catch (StopRunProcess &e) {
            this->stop(Sniper::StopRun::Peacefully);
        }
        // If return value is -1, stop the TopTask after execution
        if (-1 == this_ret) ret = -1;
    }
    if (-1 == ret) {
        return this->stop(Sniper::StopRun::Peacefully);
    }
    return true;
}
